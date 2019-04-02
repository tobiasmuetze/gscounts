#' get_covar
#' @description Get canonical form covariance matrix
#' @param timing numeric vector; 0 < \code{timing[1]} < ... < \code{timing[k]}
#' @keywords internal
get_covar <- function(timing) {
  k <- length(timing)
  # Calculate matrix with (i,j) = I_i / I_j
  covar <- sweep(x = matrix(timing, nrow = k, ncol = k), MARGIN = 2, STATS = timing, FUN = '/')
  # Make matrix symmetric
  covar[lower.tri(covar)] <- t(covar)[lower.tri(covar)]
  covar <- sqrt(covar)
  covar
}

#' get_recruittimes
#' @description Calculate the times subjects enter the study
#' @details Number of subjects in the study at study time t is given by
#' f(t)=a * t^b with  a = n / accrual_period.
#' For linear recruitment, b=1. 
#' b > 1 results is slower than linear recruitment for t < accrual_period and 
#' faster than linear recruitment for t > accrual_period. b < 1 vice versa. 
#' @param accrual_period numeric; duration of accrual period. Must be positive.
#' @param n numeric; number of subjects recruited at end of accrual period.
#' @param accrual_exponent numeric; exponent in sample size function \eqn{f(t)}. 
#' @keywords internal
get_recruittimes <- function(accrual_period, n, accrual_exponent) {
    accrual_period * (seq(0, 1, length.out = n))^(1/accrual_exponent)
}


#' samplesize_from_periods
#' @description Determine the total sample size for given accrual and study period
#' @keywords internal
samplesize_from_periods <- function(max_info, accrual_period, study_period, random_ratio,
                                    rate1, rate2, shape, accrual_speed){
  # Calculate information level depending on sample size
  body_periods <- quote({
    n1 <- ceiling(n / (1 + 1/random_ratio))
    n2 <- ceiling(n / (1 + random_ratio))
    n <- n1 + n2
    t_recruit1 <- get_recruittimes(accrual_period = accrual_period, 
                                   n = n1, 
                                   accrual_exponent = accrual_speed)
    t_recruit2 <- get_recruittimes(accrual_period = accrual_period, 
                                   n = n2, 
                                   accrual_exponent = accrual_speed)
    get_info_gsnb(rate1 = rate1, rate2 = rate2, dispersion = shape,
                  followup1 = study_period - t_recruit1,
                  followup2 = study_period - t_recruit2)
  })
  
  # Find sample size such that maximum information level is reached
  n <- uniroot(f = function(n){eval(body_periods) - max_info},
               interval = c(2*random_ratio, 1000),
               extendInt = "upX")$root
  # Update the information level with the latest sample size
  n1 <- ceiling(n / (1 + 1/random_ratio))
  n2 <- ceiling(n / (1 + random_ratio))
  n <- n1 + n2
  t_recruit1 <- get_recruittimes(accrual_period = accrual_period, 
                                 n = n1, 
                                 accrual_exponent = accrual_speed)
  t_recruit2 <- get_recruittimes(accrual_period = accrual_period, 
                                 n = n2, 
                                 accrual_exponent = accrual_speed)
  max_info <- get_info_gsnb(rate1 = rate1, rate2 = rate2, dispersion = shape,
                            followup1 = study_period - t_recruit1,
                            followup2 = study_period - t_recruit2)
  
  list(n = n, n1 = n1, n2 = n2, t_recruit1 = t_recruit1, 
       t_recruit2 = t_recruit2, max_info = max_info, 
       accrual_period = accrual_period, study_period = study_period)
}  


#' studyperiod_from_recruit
#' @description Determine the study period for given recruitment times
#' @keywords internal
studyperiod_from_recruit <- function(max_info, random_ratio, rate1, 
                                     rate2, shape, t_recruit1, t_recruit2) {
  
  if (length(t_recruit1)/length(t_recruit2) != random_ratio)
    warning(paste("k will be set to ", length(t_recruit1)/length(t_recruit2)))
  
  # Determine maximum information for given recruitment times with variable study period
  body_recruit <- quote({
    isRec1 <- (t_recruit1 <= study_period)
    isRec2 <- (t_recruit2 <= study_period)
    get_info_gsnb(rate1 = rate1, rate2 = rate2, dispersion = shape,
                  followup1 = study_period - t_recruit1[isRec1],
                  followup2 = study_period - t_recruit2[isRec2])
  })
  
  # Calculate study period for given recruitment times
  study_period <- uniroot(f = function(study_period){eval(body_recruit) - max_info},
                          interval = c(min(t_recruit1, t_recruit2), 5*max(t_recruit1, t_recruit2)),
                          extendInt = "upX")$root
  
  # Update the information level with the latest sample size
  max_info <- eval(body_recruit)
  n1 <- sum(t_recruit1 <= study_period)
  n2 <- sum(t_recruit2 <= study_period)
  n <- n1 + n2  
  list(n = n, n1 = n1, n2 = n2, t_recruit1 = t_recruit1, 
       t_recruit2 = t_recruit2, max_info = max_info, 
       study_period = study_period)
}  


#' samplesize_from_followup
#' @description Determine the sample size from a given follow-up time
#' @keywords internal
samplesize_from_followup <- function(max_info, random_ratio, rate1, 
                                     rate2, shape, followup_max) {
  
  # Determine maximum information for given fix follow-up time
  body_followmax <- quote({
    n1 <- ceiling(n / (1 + 1/random_ratio))
    n2 <- ceiling(n / (1 + random_ratio))
    n <- n1 + n2
    get_info_gsnb(rate1 = rate1, rate2 = rate2, dispersion = shape,
                  followup1 = rep(followup_max, times = n1),
                  followup2 = rep(followup_max, times = n2))
  })
  
  
  # Calculate the sample size for a fixed individual follow-up
  n <- uniroot(f = function(n){eval(body_followmax) - max_info},
               interval = c(2*random_ratio, 1000), extendInt = "upX")$root
  
  # Update the information level with the latest sample size
  n1 <- ceiling(n / (1 + 1/random_ratio))
  n2 <- ceiling(n / (1 + random_ratio))
  n <- n1 + n2
  max_info <- get_info_gsnb(rate1 = rate1, rate2 = rate2, dispersion = shape,
                            followup1 = rep(followup_max, times = n1),
                            followup2 = rep(followup_max, times = n2))
  
  
  list(n = n, n1 = n1, n2 = n2, max_info = max_info, 
       followup_max = followup_max)
}  