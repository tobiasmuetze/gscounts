#' @name get_info_gsnb
#' @title Information level for log rate ratio
#' @description Calculates the information level for the log rate ratio of
#' the negative binomial model
#' @param rate1 numeric; rate in treatment group 1
#' @param rate2 numeric; rate in treatment group 2
#' @param shape numeric; shape parameter
#' @param followup1 numeric vector; individual follow-up times in treatment group 1
#' @param followup2 numeric vector; individual follow-up times in treatment group 2
#' @return numeric; information level
#' @export
get_info_gsnb <- function(rate1, rate2, shape, followup1, followup2) {

  if (any(c(length(rate1), length(rate2), length(shape)) != 1)) stop("error in calc_info")

  reciprocal_info1 <- 1 / sum(followup1 * rate1 / (1 + shape * followup1 * rate1))
  reciprocal_info2 <- 1 / sum(followup2 * rate2 / (1 + shape * followup2 * rate2))
  1 / (reciprocal_info1 + reciprocal_info2)
}


#' @name get_timepoints_gsnb
#' @title Time point of looks
#' @description Calculate the time points of looks given the information time
#' @param rate1 numeric; rate in treatment group 1
#' @param rate2 numeric; rate in treatment group 2
#' @param shape numeric; shape parameter
#' @param t_recruit1 numeric vector; recruit (i.e. study entry) times in group 1
#' @param t_recruit2 numeric vector; recruit (i.e. study entry) times in group 2
#' @param timing numeric vector with entries in (0,1]; information times of looks
#' @param followup1 numeric vector; final individual follow-up times in treatment group 1
#' @param followup2 numeric vector; final individual follow-up times in treatment group 2
#' @return numeric; vector with the time points of the looks
get_timepoints_gsnb <- function(rate1, rate2, shape, t_recruit1, t_recruit2,
                                timing, followup1, followup2) {

  max_info <- get_info_gsnb(rate1 = rate1, rate2 = rate2, shape = shape,
                             followup1 = followup1, followup2 = followup2)

  # Information level per stage
  stage_info <- max_info * timing

  # Information at time t
  info_at_t <- quote({
    time_in_study_1 <- pmin(pmax(t - t_recruit1, 0), followup1)
    time_in_study_2 <- pmin(pmax(t - t_recruit2, 0), followup2)
    calc_info_gsnb(rate1 = rate1, rate2 = rate2, shape = shape,
                   followup1 = time_in_study_1,
                   followup2 = time_in_study_2)
  })

  # Time points at which (interim) analyses are conducted
  analysis_times <- sapply(X = stage_info,
                           FUN = function(y) {
                             uniroot(f = function(t){eval(info_at_t) - y},
                                     interval = c(0, max(followup1, followup2)),
                                     tol = .Machine$double.eps^0.75)$root
                           })
  analysis_times
}
