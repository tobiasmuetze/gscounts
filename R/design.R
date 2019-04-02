#' @name add_stopping_prob
#' @title Calculate stopping probabilities
#' @description Calculate analyses specific stopping probabilities and
#' expected information level
#' @param x Object of class gsnb
#' @return Object of class gsnb
#' @import stats
#' @keywords internal
add_stopping_prob <- function(x) {
  
  k <- length(x$timing)
  ratio_H1 <- x$rate1 / x$rate2
  log_effect <- log(ratio_H1) - log(x$ratio_H0)
  
  # Calculation of stopping probabilities and expected information levels in 3 parts:
  # (1) Probabilities for stopping for efficacy
  # (2) Probabilities for stopping for futility if futility stopping is included
  # (3) Expected information levels
  
  # (1) Probabilities for stopping for efficacy/rejecting the null 
  # hypothesis under either H0 or H1
  eff_stop <- matrix(nrow = 2, ncol = k) # first row for H0, second row for H1
  eff_stop[1, 1] <- pnorm(x$efficacy$critical[1], mean = log(x$ratio_H0) * sqrt(x$timing[1] * x$max_info))
  eff_stop[2, 1] <- pnorm(x$efficacy$critical[1], mean = log_effect * sqrt(x$timing[1] * x$max_info))
  if (is.null(x$futility)) 
    fut_critical <- rep(Inf, times = k)
  else
    fut_critical <- x$futility$critical
  ## Probs for stopping for efficacy
  for(i in 2:k) {
    eff_stop[1, i] <- cpp_pmultinorm(r = 23, 
                                     lower = c(x$efficacy$critical[1:(i-1)], -Inf), 
                                     upper = c(fut_critical[1:(i-1)], x$efficacy$critical[i]), 
                                     information = x$max_info * x$timing[1:i], 
                                     theta = log(x$ratio_H0))
    eff_stop[2, i] <- cpp_pmultinorm(r = 23, 
                                     lower = c(x$efficacy$critical[1:(i-1)], -Inf), 
                                     upper = c(fut_critical[1:(i-1)], x$efficacy$critical[i]), 
                                     information = x$max_info * x$timing[1:i], 
                                     theta = log_effect)
  }
  ## Initialize results as part of x
  x$stop_prob$efficacy <- as.data.frame(cbind(c(x$ratio_H0, ratio_H1), eff_stop, 
                                              rowSums(eff_stop)))
  names(x$stop_prob$efficacy) <- c("Rate ratio", paste0("Look ", 1:k), "Total")
  
  # (2) Probabilities for stopping for (nonbinding) futility/accepting the null
  # hypothesis under either H0 or H1
  if (!is.null(x$futility)) {
    fut_stop <- matrix(nrow = 2, ncol = k)
    fut_stop[1, 1] <- 1 - pnorm(x$futility$critical[1], mean = log(x$ratio_H0) * sqrt(x$timing[1] * x$max_info))
    fut_stop[2, 1] <- 1 - pnorm(x$futility$critical[1], mean = log_effect * sqrt(x$timing[1] * x$max_info))
    ## Probs for stopping for efficacy
    for(i in 2:k) {
      fut_stop[1, i] <- cpp_pmultinorm(r = 23, 
                                       lower = c(x$efficacy$critical[1:(i-1)], x$futility$critical[i]), 
                                       upper = c(x$futility$critical[1:(i-1)], Inf), 
                                       information = x$max_info * x$timing[1:i], 
                                       theta = log(x$ratio_H0))
      fut_stop[2, i] <- cpp_pmultinorm(r = 23, 
                                       lower = c(x$efficacy$critical[1:(i-1)], x$futility$critical[i]), 
                                       upper = c(x$futility$critical[1:(i-1)], Inf), 
                                       information = x$max_info * x$timing[1:i], 
                                       theta = log_effect)
    }
    ## Initialize results as part of x
    x$stop_prob$futility <- as.data.frame(cbind(c(x$ratio_H0, ratio_H1), fut_stop, 
                                                rowSums(fut_stop)))
    names(x$stop_prob$futility) <- c("Rate ratio", paste0("Look ", 1:k), "Total")
  }
  
  # (3) Expected information level
  if (is.null(x$futility)) 
    expected_info <- (eff_stop %*% x$timing) * x$max_info + (1 - rowSums(eff_stop)) * x$timing[k] * x$max_info
  else 
    expected_info <- ((eff_stop + fut_stop) %*% x$timing) * x$max_info + 
    (1 - rowSums(eff_stop) - rowSums(fut_stop)) * x$timing[k] * x$max_info
  
  x$expected_info <- as.data.frame(cbind(c(x$ratio_H0, ratio_H1), expected_info))
  names(x$expected_info) <- c("Rate ratio", "E[I]")
  
  # Return object of class gsnb
  x
}






#' @name design_gsnb
#' @title Group sequential design with negative binomial outcomes
#' @description Design a group sequential trial with negative binomial outcomes
#' @details 
#' Denote  \eqn{\mu_1} and \eqn{\mu_2} the event rates in treatment groups 1 and 2.
#' This function considers smaller event rates to be better. 
#' The statistical hypothesis testing problem of interest is 
#' \deqn{H_0: \frac{\mu_1}{\mu_2} \ge \delta  vs.  H_1: \frac{\mu_1}{\mu_2} < \delta,}
#' with \eqn{\delta=}\code{ratio_H0}.
#' Non-inferiority of treatment group 1 compared to treatment group 2 is tested for \eqn{\delta\in (1,\infty)}.
#' Superiority of treatment group 1 over treatment group 2 is tested for \eqn{\delta \in (0,1]}.
#' The calculation of the efficacy and (non-)binding futility boundaries are performed
#' under the hypothesis \eqn{H_0: \frac{\mu_1}{\mu_2}= \delta} and 
#' under the alternative \eqn{H_1: \frac{\mu_1}{\mu_2} = }\code{rate1} / \code{rate2}.
#' 
#' The argument `accrual_speed` is used to adjust the accrual speed.
#' Number of subjects in the study at study time t is given by
#' \eqn{f(t)=a * t^b} with  \eqn{a = n / accrual_period} and \eqn{b=accrual_speed}  
#' For linear recruitment, \eqn{b=1}. 
#' \eqn{b > 1} results is slower than linear recruitment for \eqn{t < accrual_period} and 
#' faster than linear recruitment for \eqn{t > accrual_period}. Vice verse for \eqn{b < 1}. 
#' @param rate1 numeric; assumed rate of treatment group 1 in the alternative
#' @param rate2 numeric; assumed rate of treatment group 2 in the alternative
#' @param dispersion numeric; dispersion (shape) parameter of negative binomial distribution
#' @param power numeric; target power of group sequential design
#' @param timing numeric vector; 0 < \code{timing[1]} < ... < \code{timing[K]} = 1
#' with \code{K} the number of analyses, i.e. (K-1) interim analyses and final analysis.
#' When the timing of efficacy and futility analyses differ, timing should not be defined.
#' Instead, the arguments \code{timing_eff} and \code{timing_fut} have to be used to specify
#' the timing of the efficacy and futility analyses, respectively. 
#' @param esf function; error spending function
#' @param ratio_H0 numeric; positive number denoting the rate ratio \eqn{\mu_1/\mu_2}
#' under the null hypothesis, i.e. the non-inferiority or superiority margin
#' @param sig_level numeric; Type I error / significance level
#' @param random_ratio numeric; randomization ratio n1/n2
#' @param futility character; either \code{"binding"}, \code{"nonbinding"}, or 
#' \code{NULL} for binding, nonbinding, or no futility boundaries
#' @param esf_futility function; futility error spending function
#' @param t_recruit1 numeric vector; recruit (i.e. study entry) times in group 1
#' @param t_recruit2 numeric vector; recruit (i.e. study entry) times in group 2
#' @param study_period numeric; study duration;
#' to be set when follow-up times are not identical between subjects, NULL otherwise
#' @param accrual_period numeric; accrual period
#' @param followup_max numeric; maximum exposure time of a subject; 
#' to be set when follow-up times are to be equal for each subject, NULL otherwise
#' @param accrual_speed numeric; determines accrual speed; values larger than 1
#' result in accrual slower than linear; values between 0 and 1 result in accrual 
#' faster than linear.  
#' @param ... further arguments. Will be passed to the error spending function.
#' @return A list with class "gsnb" containing the following components:
#' \item{rate1}{as input}
#' \item{rate2}{as input}
#' \item{dispersion}{as input}
#' \item{power}{as input}
#' \item{timing}{as input}
#' \item{ratio_H0}{as input}
#' \item{ratio_H1}{ratio \code{rate1}/\code{rate2}}
#' \item{sig_level}{as input}
#' \item{random_ratio}{as input}
#' \item{power_fix}{power of fixed design}
#' \item{expected_info}{list; expected information under \code{ratio_H0} and \code{ratio_H1}}
#' \item{efficacy}{list; contains the elements \code{esf} (type I error spending function), 
#'                \code{spend} (type I error spend at each look), and 
#'                \code{critical} (critical value for efficacy testing)}
#' \item{futility}{list; only part of the output if argument \code{futility} is
#'                 defined in the input. Contains the elements \code{futility} 
#'                 (input argument \code{futility}), \code{esf} 
#'                 (type II error spending function), \code{spend} (type II 
#'                 error spend at each look), and \code{critical} (critical 
#'                 value for futility testing)}
#' \item{stop_prob}{list; contains the element \code{efficacy} with the probabilities
#'                  for stopping for efficacy and, if futility bounds are calculated, 
#'                  the element \code{futility} with the probabilities for stopping 
#'                  for futility}               
#' \item{t_recruit1}{as input}
#' \item{t_recruit2}{as input}
#' \item{study_period}{as input}
#' \item{followup_max}{as input}
#' \item{max_info}{maximum information}
#' \item{calendar}{calendar times of data looks; only calculated when exposure times are not identical}
#' @references Mütze, T., Glimm, E., Schmidli, H., & Friede, T. (2018). 
#' Group sequential designs for negative binomial outcomes. 
#' Statistical Methods in Medical Research, \url{https://doi.org/10.1177/0962280218773115}.
#' @examples
#' # Calculate the sample sizes for a given accrual period and study period (without futility)
#' out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, dispersion = 5, 
#'                    power = 0.8, timing = c(0.5, 1), esf = obrien,
#'                    ratio_H0 = 1, sig_level = 0.025,
#'                    study_period = 3.5, accrual_period = 1.25, random_ratio = 1)
#' out
#' 
#' # Calculate the sample sizes for a given accrual period and study period with binding futility
#' out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, dispersion = 5, 
#'                    power = 0.8, timing = c(0.5, 1), esf = obrien,
#'                    ratio_H0 = 1, sig_level = 0.025, study_period = 3.5, 
#'                    accrual_period = 1.25, random_ratio = 1, futility = "binding", 
#'                    esf_futility = obrien)
#' out
#' 
#' 
#' # Calculate study period for given recruitment times
#' expose <- seq(0, 1.25, length.out = 1042)
#' out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, dispersion = 5, 
#'                    power = 0.8, timing = c(0.5, 1), esf = obrien,
#'                    ratio_H0 = 1, sig_level = 0.025, t_recruit1 = expose, 
#'                    t_recruit2 = expose, random_ratio = 1)
#' out
#'
#' # Calculate sample size for a fixed exposure time
#' out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, dispersion = 5, 
#'                    power = 0.8, timing = c(0.5, 1), esf = obrien,
#'                    ratio_H0 = 1, sig_level = 0.025,
#'                    followup_max = 0.5, random_ratio = 1)
#'                    
#' # Different timing for efficacy and futility analyses
#'  design_gsnb(rate1 = 1, rate2 = 2, dispersion = 5,
#'              power = 0.8, esf = obrien,
#'              ratio_H0 = 1, sig_level = 0.025, study_period = 3.5,
#'              accrual_period = 1.25, random_ratio = 1, futility = "binding",
#'              esf_futility = pocock, 
#'              timing_eff = c(0.8, 1),
#'              timing_fut = c(0.2, 0.5, 1))                    
#' @export
design_gsnb <- function(rate1, rate2, dispersion, ratio_H0 = 1, random_ratio = 1, 
                        power, sig_level, timing, esf = obrien,
                        esf_futility = NULL, futility = NULL,
                        t_recruit1 = NULL, t_recruit2 = NULL, study_period = NULL, 
                        accrual_period = NULL, followup_max = NULL, 
                        accrual_speed = 1, ...) {
  
  # Read out separate input of timing for efficacy and futility stopping
  timing_eff <- list(...)$timing_eff
  timing_fut <- list(...)$timing_fut

  # Error check for timing variables  
  if (missing(timing) & !is.null(timing_eff) & !is.null(timing_fut)) {
    timing <- sort(unique(c(timing_eff, timing_fut)))
  } else if (missing(timing) & (is.null(timing_eff) | is.null(timing_fut))) {
    stop("Argument 'timing' is missing and 'timing_eff' or 'timing_fut' is missing, too.")
  } else {
    if (!is.null(timing_eff) | !is.null(timing_fut)) {
      warning("Argument 'timing' is specified; 'timing_eff' and 'timing_fut' will be overwritten.")
    }
    timing_eff <- timing_fut <- timing
  }
    
  isNull_fm <- is.null(followup_max)
  isNull_sp <- is.null(study_period)
  isNull_ap <- is.null(accrual_period)
  isNull_t1 <- is.null(t_recruit1)
  isNull_t2 <- is.null(t_recruit2)
  
  # Error check for accrual speed
  if (accrual_speed <= 0) stop("accrual_speed must be positive")
  
  # Initialize object of class gsnb
  x <- list(rate1 = rate1, rate2 = rate2, dispersion = dispersion, 
            ratio_H0 = ratio_H0, power = power, sig_level = sig_level, 
            timing = timing, random_ratio = random_ratio)
  # Efficacy spending
  esf_eff_out <- esf(t = timing_eff, 
                     sig_level = sig_level, ...)
  spend_eff <- c(esf_eff_out[1], diff(esf_eff_out))
  # Futility spending
  if (!is.null(futility)) {
    esf_fut_out <- esf_futility(t = timing_fut, sig_level = 1 - power, ...)
    spend_fut <- c(esf_fut_out[1], diff(esf_fut_out))
    df_fut <- data.frame(timing = timing_fut, spend_fut = spend_fut)
    
    df_eff <- data.frame(timing = timing_eff, spend_eff = spend_eff)
    df_merge <- merge(df_eff, df_fut, all = TRUE)
    df_merge[is.na(df_merge)] <- 0
    
    x$futility <- list(esf = esf_futility,
                       spend = df_merge$spend_fut,
                       type = futility)
    x$efficacy <- list(esf = esf, 
                       spend = df_merge$spend_eff)
  } else {
    x$efficacy <- list(esf = esf, 
                       spend = spend_eff)
  }
  x$ratio_H1 <- rate1 / rate2
  class(x) <- "gsnb"
  # Perform validity checks
  check_gsnb(x)
  
  # Add maximum information to x
  x <- add_maxinfo(x)
  
  # Calcualte the missing design characteristics
  # 1. Sample size in the case of a fixed followup.
  # 2. Study period for given recruitment times
  # 3. Sample size for fixed accrual and study period (using uniform recruitment)
  if (all(!isNull_fm, isNull_sp, isNull_t1, isNull_t2)) {
    # 1. Calculate the sample size for a fixed individual follow-up
    out <- samplesize_from_followup(max_info = x$max_info, random_ratio = random_ratio,
                                    rate1 = rate1, rate2 = rate2, shape = dispersion,
                                    followup_max = followup_max)
    # Add calendar time of data look if accrual period or recruitment times are defined
    if (!isNull_ap) {
      out$t_recruit1 <- get_recruittimes(accrual_period = accrual_period, 
                                         n = out$n1, 
                                         accrual_exponent = accrual_speed)
      out$t_recruit2 <- get_recruittimes(accrual_period = accrual_period, 
                                         n = out$n2, 
                                         accrual_exponent = accrual_speed)
      
      x$calendar <- get_calendartime_gsnb(rate1 = rate1, rate2 = rate2, dispersion = dispersion, 
                                          t_recruit1 = out$t_recruit1, t_recruit2 = out$t_recruit2,
                                          timing = timing, 
                                          followup1 = rep(followup_max, times = out$n1),
                                          followup2 = rep(followup_max, times = out$n2))
    }
  } else if (all(isNull_sp, isNull_fm, isNull_ap, !isNull_t1, !isNull_t1)) {
    # 2. Calculate study period for given recruitment times
    out <- studyperiod_from_recruit(max_info = x$max_info, random_ratio = random_ratio,
                                    rate1 = rate1, rate2 = rate2, shape = dispersion,
                                    t_recruit1 = t_recruit1, t_recruit2 = t_recruit2)
    x$calendar <- get_calendartime_gsnb(rate1 = rate1, rate2 = rate2, dispersion = dispersion, 
                                        t_recruit1 = out$t_recruit1, t_recruit2 = out$t_recruit2,
                                        timing = timing, followup1 = out$study_period - out$t_recruit1,
                                        followup2 = out$study_period - out$t_recruit2)
  } else if (all(!isNull_sp, !isNull_ap, isNull_fm, isNull_t1, isNull_t1)) {
    # 3. Calculate sample size for fixed accrual and study period (using uniform recruitment)
    out <- samplesize_from_periods(max_info = x$max_info, accrual_period = accrual_period,
                                   study_period = study_period,
                                   random_ratio = random_ratio,
                                   rate1 = rate1, rate2 = rate2, shape = dispersion,
                                   accrual_speed = accrual_speed)
    x$calendar <- get_calendartime_gsnb(rate1 = rate1, rate2 = rate2, dispersion = dispersion, 
                                        t_recruit1 = out$t_recruit1, 
                                        t_recruit2 = out$t_recruit2,
                                        timing = timing, 
                                        followup1 = study_period - out$t_recruit1,
                                        followup2 = study_period - out$t_recruit2)
  } else {
    stop("No appropriate combination of input arguments is defined")
  }
  
  
  # Remove maxium information from x since it is recalculated in out
  x$max_info <- NULL
  x <- c(x, out)
  class(x) <- "gsnb"
  # Add stopping probabilities to x
  x <- add_stopping_prob(x)
  # Add power of fixed design to x
  log_effect <- log(x$rate1 / x$rate2) - log(x$ratio_H0)
  x$power_fix <- pnorm(qnorm(sig_level), mean = log_effect * sqrt(x$max_info))
  
  x
}

