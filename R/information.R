#' @name get_info_gsnb
#' @title Information level for log rate ratio
#' @description Calculates the information level for the log rate ratio of
#' the negative binomial model
#' @param rate1 numeric; rate in treatment group 1
#' @param rate2 numeric; rate in treatment group 2
#' @param dispersion numeric; dispersion (shape) parameter of negative binomial distribution
#' @param followup1 numeric vector; individual follow-up times in treatment group 1
#' @param followup2 numeric vector; individual follow-up times in treatment group 2
#' @return numeric; information level
#' @examples 
#' # Calculates information level for case of 10 subjects per group
#' # Follow-up times of subjects in each group range from 1 to 3
#' get_info_gsnb(rate1 = 0.1,
#'               rate2 = 0.125,
#'               dispersion = 4, 
#'               followup1 = seq(1, 3, length.out = 10), 
#'               followup2 = seq(1, 3, length.out = 10))
#' @export
get_info_gsnb <- function(rate1, rate2, dispersion, followup1, followup2) {

  if (any(c(length(rate1), length(rate2), length(dispersion)) != 1)) stop("error in calc_info")

  reciprocal_info1 <- 1 / sum(followup1 * rate1 / (1 + dispersion * followup1 * rate1))
  reciprocal_info2 <- 1 / sum(followup2 * rate2 / (1 + dispersion * followup2 * rate2))
  1 / (reciprocal_info1 + reciprocal_info2)
}


#' @name get_calendartime_gsnb
#' @title Calendar time of data looks
#' @description Calculate the calendar time of looks given the information time
#' @param rate1 numeric; rate in treatment group 1
#' @param rate2 numeric; rate in treatment group 2
#' @param dispersion numeric; dispersion (shape) parameter of negative binomial distribution
#' @param t_recruit1 numeric vector; recruit (i.e. study entry) times in group 1
#' @param t_recruit2 numeric vector; recruit (i.e. study entry) times in group 2
#' @param timing numeric vector with entries in (0,1]; information times of data looks
#' @param followup1 numeric vector; final individual follow-up times in treatment group 1
#' @param followup2 numeric vector; final individual follow-up times in treatment group 2
#' @return numeric; vector with calendar time of data looks
#' @import stats
#' @examples 
#' # Calendar time at which 50%, 75%, and 100% of the maximum information is attained
#' # 100 subjects in each group are recruited uniformly over 1.5 years
#' # Study ends after two years, i.e. follow-up times vary between 2 and 0.5 years 
#' get_calendartime_gsnb(rate1 = 0.1, 
#'                       rate2 = 0.125, 
#'                       dispersion = 5, 
#'                       t_recruit1 = seq(0, 1.5, length.out = 100), 
#'                       t_recruit2 = seq(0, 1.5, length.out = 100),
#'                       timing = c(0.5, 0.75, 1),
#'                       followup1 = seq(2, 0.5, length.out = 100),
#'                       followup2 = seq(2, 0.5, length.out = 100)) 
#' @export
get_calendartime_gsnb <- function(rate1, rate2, dispersion, t_recruit1, t_recruit2,
                                timing, followup1, followup2) {

  max_info <- get_info_gsnb(rate1 = rate1, rate2 = rate2, dispersion = dispersion,
                             followup1 = followup1, followup2 = followup2)

  # Information level per stage
  stage_info <- max_info * timing

  # Information at time t
  info_at_t <- quote({
    time_in_study_1 <- pmin(pmax(t - t_recruit1, 0), followup1)
    time_in_study_2 <- pmin(pmax(t - t_recruit2, 0), followup2)
    get_info_gsnb(rate1 = rate1, rate2 = rate2, 
                  dispersion = dispersion,
                   followup1 = time_in_study_1,
                   followup2 = time_in_study_2)
  })

  # Calendar times at which (interim) analyses are conducted
  analysis_times <- sapply(X = stage_info,
                           FUN = function(y) {
                             uniroot(f = function(t){eval(info_at_t) - y},
                                     interval = c(0.001, max(t_recruit1 + followup1, t_recruit2 + followup2)),
                                     tol = .Machine$double.eps^0.75)$root
                           })
  analysis_times
}
