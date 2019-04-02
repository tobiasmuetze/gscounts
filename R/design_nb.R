#' @name design_nb
#' @title Clinical trials with negative binomial outcomes
#' @description Design a clinical trial with negative binomial outcomes
#' @param rate1 numeric; assumed rate of treatment group 1 in the alternative
#' @param rate2 numeric; assumed rate of treatment group 2 in the alternative
#' @param dispersion numeric; dispersion (shape) parameter of negative binomial distribution
#' @param power numeric; target power
#' @param ratio_H0 numeric; positive number denoting the rate ratio rate_1/rate_2
#' under the null hypothesis, i.e. the non-inferiority or superiority margin
#' @param sig_level numeric; Type I error / significance level
#' @param random_ratio numeric; randomization ratio n1/n2
#' @param t_recruit1 numeric vector; recruit (i.e. study entry) times in group 1
#' @param t_recruit2 numeric vector; recruit (i.e. study entry) times in group 2
#' @param study_period numeric; study duration
#' @param accrual_period numeric; accrual period
#' @param accrual_speed numeric; determines accrual speed; values larger than 1
#' result in accrual slower than linear; values between 0 and 1 result in accrual 
#' faster than linear.  
#' @param followup_max numeric; maximum exposure time of a patient
#' @return A list containing the following components:
#' \item{rate1}{as input}
#' \item{rate2}{as input}
#' \item{dispersion}{as input}
#' \item{power}{as input}
#' \item{ratio_H0}{as input}
#' \item{ratio_H1}{ratio \code{rate1}/\code{rate2}}
#' \item{sig_level}{as input}
#' \item{random_ratio}{as input}
#' \item{t_recruit1}{as input}
#' \item{t_recruit2}{as input}
#' \item{study_period}{as input}
#' \item{followup_max}{as input}
#' \item{max_info}{maximum information}
#' @examples
#' # Calculate sample size for given accrual period and study duration assuming uniformal accrual
#' out <- design_nb(rate1 = 0.0875, rate2 = 0.125, dispersion = 5, power = 0.8,
#'                  ratio_H0 = 1, sig_level = 0.025,
#'                  study_period = 4, accrual_period = 1, random_ratio = 2)
#' out
#' 
#' # Calculate sample size for a fixed exposure time of 0.5 years
#' out <- design_nb(rate1 = 4.2, rate2 = 8.4, dispersion = 3, power = 0.8,
#'                  ratio_H0 = 1, sig_level = 0.025,
#'                  followup_max = 0.5, random_ratio = 2)
#' out
#' 
#' # Calculate study period for given recruitment time
#' t_recruit1 <- seq(0, 1.25, length.out = 1200)
#' t_recruit2 <- seq(0, 1.25, length.out = 800)
#' out <- design_nb(rate1 = 0.0875, rate2 = 0.125, dispersion = 5, power = 0.8,
#'                  ratio_H0 = 1, sig_level = 0.025,
#'                  t_recruit1 = t_recruit1, t_recruit2 = t_recruit2)
#' @import stats
#' @export
design_nb <- function(rate1, rate2, dispersion, power, ratio_H0 = 1, sig_level,
                      random_ratio = 1, t_recruit1 = NULL,
                      t_recruit2 = NULL, study_period = NULL, accrual_period = NULL,
                      followup_max = NULL, 
                      accrual_speed = 1) {
  
  arguments <- as.list(environment())
  
  # Error check for accrual speed
  if (accrual_speed <= 0) stop("accrual_speed must be positive")
  
  # Calculate maximum information required to obtain power
  max_info <- (qnorm(1-sig_level) + qnorm(power))^2 / log(rate1 / rate2 / ratio_H0)^2
  
  isNull_fm <- is.null(followup_max)
  isNull_sp <- is.null(study_period)
  isNull_ap <- is.null(accrual_period)
  isNull_t1 <- is.null(t_recruit1)
  isNull_t2 <- is.null(t_recruit2)
  
  # Calcualte the missing design value
  if (all(!isNull_fm, isNull_sp, isNull_ap)) {
    # Calculate the sample size for a fixed individual follow-up
    out <- samplesize_from_followup(max_info = max_info, random_ratio = random_ratio, 
                                    rate1 = rate1, rate2 = rate2, shape = dispersion, 
                                    followup_max = followup_max)
  } else if (all(isNull_sp, isNull_fm, isNull_ap, !isNull_t1, !isNull_t1)) {
    # Calculate study period for given recruitment times
    out <- studyperiod_from_recruit(max_info = max_info, random_ratio = random_ratio, 
                                    rate1 = rate1, rate2 = rate2, shape = dispersion, 
                                    t_recruit1 = t_recruit1, t_recruit2 = t_recruit2)
  } else if (all(!isNull_sp, !isNull_ap, isNull_fm, isNull_t1, isNull_t1)) {
    # Calculate sample size for fixed accrual and study period (using uniform recruitment)
    out <- samplesize_from_periods(max_info = max_info, accrual_period = accrual_period, 
                                   study_period = study_period, 
                                   random_ratio = random_ratio,
                                   rate1 = rate1, rate2 = rate2, shape = dispersion,
                                   accrual_speed = accrual_speed)
  } else {
    stop("No appropriate combination of input arguments is defined")
  }
  
  
  out <- c(arguments[c("rate1", "rate2", "dispersion", "power", "ratio_H0", "sig_level")], out)
  class(out) <- "nb"
  out
}
