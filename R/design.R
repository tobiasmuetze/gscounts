#' @name get_rejectprob_gsnb
#' @title Calcualte rejection probabilities
#' @description Calcualte analyses specific rejection probabilities and
#' expected information level, study duration, and sample size
#' @param rate_ratio The rate ratio for which the rejection probability is calculated.
#' @param ratio_H0 The rate ratio in the null hypothesis against which the test is performed.
#' @param critical The vector of critical values used for the test. Must have
#' the same length as argument \code{timing}.
#' @param max_info The maximum information of the study.
#' @param timing The vector of information times at which the tests are performed.
#' @param n_vec An optional vector with \code{n_vec[i]} indicating the number of
#' subjects recruited at information time \code{timing[i]}.
#' Must have the same length as argument \code{timing}.
#' @param study_period_vec An optional vector with \code{study_period_vec[i]}
#' indicating the study time at information time \code{timing[i]}.
#' Must have the same length as argument \code{timing}.
#' @return data frame with rejection (i.e. boundary crossing) probabilities
#' @import mvtnorm
#' @import utils
#' @keywords internal
get_rejectprob_gsnb <- function(rate_ratio, ratio_H0 = 1, critical, max_info, timing,
                                n_vec = NULL, study_period_vec = NULL) {

  k <- length(timing)
  effect_size <- log(rate_ratio) - log(ratio_H0)
  covar <- get_covar(timing)

  # Calculate analysis specific rejection probabilities
  reject_prob <- numeric(k)
  reject_prob[1] <- pnorm(critical[1], mean = effect_size * sqrt(timing[1] * max_info))
  for(i in 2:k) {
    lower <- c(critical[1:(i-1)], -Inf)
    upper <- c(rep(Inf, times = i-1), critical[i])
    mean_vec <- effect_size * sqrt(timing[1:i] * max_info)
    reject_prob[i] <- pmvnorm(lower = lower, upper = upper,
            mean = mean_vec, sigma = covar[1:i, 1:i])[1]
  }

  # Expected information level
  expected_info <- sum(reject_prob * timing * max_info) + (1 - sum(reject_prob)) * max_info
  # Initialize output vector
  out <- c(rate_ratio, reject_prob, sum(reject_prob), expected_info)
  names(out) <- c("Rate ratio", paste0("Analysis ", 1:k), "Total", "E[I]")

  # Calculate expected sample size and study duration
  if (!is.null(n_vec)) {
    expected_ss <- sum(reject_prob * n_vec) + (1 - sum(reject_prob)) * tail(n_vec, n = 1)
    out <- c(out, expected_ss)
    names(out)[length(out)] <- "E[n]"
  }
  if (!is.null(study_period_vec)) {
    expected_study_period <- sum(reject_prob * study_period_vec) + (1 - sum(reject_prob)) * tail(study_period_vec, n = 1)
    out <- c(out, expected_study_period)
    names(out)[length(out)] <- "E[t]"
  }
  # Return value
  out
}


#' @name design_gsnb
#' @title Design a group sequential trial with negative binomial data
#' @description Design a group sequential trial with negative binomial data
#' @param rate1 numeric; assumed rate of treatment group 1 in the alternative
#' @param rate2 numeric; assumed rate of treatment group 2 in the alternative
#' @param shape numeric; shape parameter
#' @param power_gs numeric; target power of group sequential design
#' @param timing numeric vector; 0 < \code{timing[1]} < ... < \code{timing[K]} = 1
#' with \code{K} the number of analyses, i.e. (K-1) interim analyses and final analysis
#' @param esf function; error spending function
#' @param ratio_H0 numeric; postive number denoting the rate ratio rate_1/rate_2
#' under the null hypothesis, i.e. the non-inferiority or superiority margin
#' @param sig_level numeric; Type I error / significance level
#' @param random_ratio numeric; randomization ratio n1/n2
#' @param t_recruit1 numeric vector; recruit (i.e. study entry) times in group 1
#' @param t_recruit2 numeric vector; recruit (i.e. study entry) times in group 2
#' @param study_period numeric; study duration
#' @param accrual_period numeric; accrual period
#' @param followup_max numeric; maximum exposure time of a patient
#' @return A list with class "GSNBdesign" containing the following components:
#' \item{rate1}{as input}
#' \item{rate2}{as input}
#' \item{shape}{as input}
#' \item{power_gs}{as input}
#' \item{timing}{as input}
#' \item{ratio_H0}{as input}
#' \item{ratio_H1}{ratio \code{rate1}/\code{rate2}}
#' \item{sig_level}{as input}
#' \item{random_ratio}{as input}
#' \item{t_recruit1}{as input}
#' \item{t_recruit2}{as input}
#' \item{study_period}{as input}
#' \item{followup_max}{as input}
#' \item{max_info}{ maximum information}
#' \item{power_fix}{power of fixed design}
#' \item{critical}{critical values for each look}
#' \item{reject_prob}{data frame with rejection (i.e. boundary crossing) probabilities}
#' @examples
#' # Calculate the sample sizes for a given accrual period and study period
#' out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
#'                    timing = c(0.5, 1), esf = esf_obrien,
#'                    ratio_H0 = 1, sig_level = 0.025,
#'                    study_period = 3.5, accrual_period = 1.25, random_ratio = 1)
#' out
#'
#' # Calculate study period for given recruitment times
#' expose <- seq(0, 1.25, length.out = 1042)
#' out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
#'                    timing = c(0.5, 1), esf = esf_obrien,
#'                    ratio_H0 = 1, sig_level = 0.025,
#'                    t_recruit1 = expose, t_recruit2 = expose, random_ratio = 1)
#' out
#'
#' # Calculate sample size for a fixed exposure time
#' out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
#'                    timing = c(0.5, 1), esf = esf_obrien,
#'                    ratio_H0 = 1, sig_level = 0.025,
#'                    followup_max = 0.5, random_ratio = 1)
#' @export
design_gsnb <- function(rate1, rate2, shape, power_gs, timing, esf = esf_obrien,
                        ratio_H0 = 1, sig_level, random_ratio = 1, t_recruit1 = NULL,
                        t_recruit2 = NULL, study_period = NULL, accrual_period = NULL,
                        followup_max = NULL) {

  arguments <- as.list(environment())
  K <- length(timing)
  # Calculate maximum information required to obtain power power_gs
  power_out <-  power_gsnb(ratio_H1 = rate1 / rate2, power_gs = power_gs, timing = timing,
                           esf = esf, ratio_H0 = ratio_H0, sig_level = sig_level)
  max_info <- power_out$max_info
  critical <- power_out$critical
  power_fix <- power_out$power_fix

  # # Determine maximum information for given fix follow-up time
  # body_followmax <- quote({
  #   n1 <- ceiling(n / (1 + 1/random_ratio))
  #   n2 <- ceiling(n / (1 + random_ratio))
  #   n <- n1 + n2
  #   get_info_gsnb(rate1 = rate1, rate2 = rate2, shape = shape,
  #                 followup1 = rep(followup_max, times = n1),
  #                 followup2 = rep(followup_max, times = n2))
  # })
  # # Determine maximum information for given recruitment times with variable study period
  # body_recruit <- quote({
  #   isRec1 <- (t_recruit1 <= study_period)
  #   isRec2 <- (t_recruit2 <= study_period)
  #   get_info_gsnb(rate1 = rate1, rate2 = rate2, shape = shape,
  #                 followup1 = study_period - t_recruit1[isRec1],
  #                 followup2 = study_period - t_recruit2[isRec2])
  # })
  # # Determine maximum information for given accural and study period
  # # with variable sample size from accural
  # body_accural <- quote({
  #   n1 <- ceiling(n / (1 + 1/random_ratio))
  #   n2 <- ceiling(n / (1 + random_ratio))
  #   n <- n1 + n2
  #   t_recruit1 <- seq(0, accrual_period, length.out = n1)
  #   t_recruit2 <- seq(0, accrual_period, length.out = n2)
  #   get_info_gsnb(rate1 = rate1, rate2 = rate2, shape = shape,
  #                 followup1 = study_period - t_recruit1,
  #                 followup2 = study_period - t_recruit2)
  # })


  isNull_fm <- is.null(followup_max)
  isNull_sp <- is.null(study_period)
  isNull_ap <- is.null(accrual_period)
  isNull_t1 <- is.null(t_recruit1)
  isNull_t2 <- is.null(t_recruit2)

  # Calcualte the missing design value
  if (all(!isNull_fm, isNull_sp, isNull_ap)) {
    # Calculate the sample size for a fixed individual follow-up
    # n <- uniroot(f = function(n){eval(body_followmax) - max_info},
    #              interval = c(2*random_ratio, 1000), extendInt = "upX")$root
    # max_info <- eval(body_followmax)
    out <- samplesize_from_followup(max_info = max_info, random_ratio = random_ratio, 
                                    rate1 = rate1, rate2 = rate2, shape = shape, 
                                    followup_max = followup_max)
  } else if (all(isNull_sp, isNull_fm, isNull_ap, !isNull_t1, !isNull_t1)) {
    # # Calculate study period for given recruitment times
    # study_period <- uniroot(f = function(study_period){eval(body_recruit) - max_info},
    #                         interval = c(min(t_recruit1, t_recruit2), 3*max(t_recruit1, t_recruit2)),
    #                         extendInt = "upX")$root
    # max_info <- eval(body_recruit)
    # n1 <- sum(t_recruit1 <= study_period)
    # n2 <- sum(t_recruit2 <= study_period)
    # n <- n1 + n2
    out <- studyperiod_from_recruit(max_info = max_info, random_ratio = random_ratio, 
                                    rate1 = rate1, rate2 = rate2, shape = shape, 
                                    t_recruit1 = t_recruit1, t_recruit2 = t_recruit2)
  } else if (all(!isNull_sp, !isNull_ap, isNull_fm, isNull_t1, isNull_t1)) {
    # # Calculate sample size for fixed accural and study period (using uniform recruitment)
    # n <- uniroot(f = function(n){eval(body_accural) - max_info},
    #              interval = c(2*random_ratio, 1000),
    #              extendInt = "upX")$root
    # max_info <- eval(body_accural)
    out <- samplesize_from_periods(max_info = max_info, accrual_period = accrual_period,
                                   study_period = study_period,
                                   random_ratio = random_ratio,
                                   rate1 = rate1, rate2 = rate2, shape = shape)
  } else {
    stop("No appropriate combination of input arguments is defined")
  }


  # Get analysis specific rejection probabilities and expected values
  reject_prob_H1 <- get_rejectprob_gsnb(rate_ratio = rate1 / rate2,
                                        ratio_H0 = ratio_H0, critical = critical,
                                        max_info = max_info, timing = timing)
  reject_prob_H0 <- get_rejectprob_gsnb(rate_ratio = ratio_H0,
                                        ratio_H0 = ratio_H0, critical = critical,
                                        max_info = max_info, timing = timing)
  reject_prob <- rbind(reject_prob_H0, reject_prob_H1)
  rownames(reject_prob) <- NULL

  # out <- c(arguments[c("rate1", "rate2", "shape", "power_gs", "timing", "ratio_H0", "sig_level", "esf")],
  #          list(power_fix = power_fix, n = n, n1 = n1, n2 = n2, max_info = max_info, critical = power_out$critical,
  #               t_recruit1 = t_recruit1, t_recruit2 = t_recruit2,
  #               study_period = study_period, accrual_period = accrual_period,
  #               reject_prob = reject_prob))
  out <- c(arguments[c("rate1", "rate2", "shape", "power_gs", "timing", "ratio_H0", "sig_level", "esf")], 
           out, list(reject_prob = reject_prob, critical = power_out$critical))
  
  class(out) <- "GSNBdesign"
  out
}

