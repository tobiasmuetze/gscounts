#' Next steps
#' Have three different functions in design_gsnb for the different method. This must be the next step.
#' Then: Output with similar content to design_gsnb
# rm(list = ls())
# #devtools::load_all()
# library(gsnbin)
# # power_gsnb
# # Calculate power for given maximum information
# out <- power_gsnb(ratio_H1 = 0.7, timing = c(0.5, 1), max_info = 69.26,
#                   sig_level = 0.025, esf = esf_pocock)
# out$critical
# print(out)
# plot(out)
# # '
# # Calculate maximum information for a given power
# out <- power_gsnb(ratio_H1 = 0.7, timing = c(0.5, 1), power = 0.8,
#                   sig_level = 0.025, esf = esf_pocock)
# out$critical
# print(out)
# plot(out)
# '
# '
#'
#'
#'
#'
#' # design_gsnb
# Determine the sample sizes for a given accrual period and study period
# out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
#                    timing = c(0.5, 1), esf = esf_obrien,
#                    ratio_H0 = 1, sig_level = 0.025,
#                    study_period = 3.5, accrual_period = 1.25, random_ratio = 1)
# out
#
# # Calculate study period for given recruitment times
# expose <- seq(0, 1.25, length.out = 1042)
# out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
#                    timing = c(0.5, 1), esf = esf_obrien,
#                    ratio_H0 = 1, sig_level = 0.025,
#                    t_recruit1 = expose, t_recruit2 = expose, random_ratio = 1)
# out
#
# Calculate sample size for a fixed exposure time
# out <- design_gsnb(rate1 = 4.2, rate2 = 8.4, shape = 3, power_gs = 0.8,
#                    timing = c(0.5, 1), esf = esf_obrien,
#                    ratio_H0 = 1, sig_level = 0.025,
#                    followup_max = 0.5, random_ratio = 1)
# out$t_recruit1
#'
#'
# get_timepoints_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5,
#                      t_recruit1 = out$t_recruit1,
#                      t_recruit2 = out$t_recruit2,
#                      followup1 = out$study_period - out$t_recruit1,
#                      followup2 = out$study_period - out$t_recruit2,
#                      timing = c(0.5, 1))
#
#' out <- design_gsnb(rate1 = 4.2, rate2 = 8.4, shape = 3, power_gs = 0.8,
#'                    timing = c(0.5, 1), esf = esf_obrien,
#'                    ratio_H0 = 1, sig_level = 0.025,
#'                    followup_max = 0.5, random_ratio = 1)
#' critical <- out$critical
#' max_info <- out$max_info
#' timing <- out$timing
#'
#'
#' power_gsnb(ratio_H1 = 0.7,
#'            timing = c(0.5, 1), max_info = 61.94334, sig_level = 0.025, esf = esf_obrien)
#'
#' design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
#'             timing = c(0.5, 1), esf = esf_obrien,
#'             ratio_H0 = 1, sig_level = 0.025,
#'             t_recruit1 = seq(0, 1.25, length.out = 1041),
#'             t_recruit2 = seq(0, 1.25, length.out = 1041),
#'             random_ratio = 1)
#'
#'
#' rate1 = 0.0875
#' rate2 = 0.125
#' shape = 5
#' power_gs = 0.8
#' timing = c(0.5, 1)
#' esf = esf_obrien
#' ratio_H0 = 1
#' sig_level = 0.025
#' t_recruit1 = seq(0, 1.25, length.out = 1043)
#' t_recruit2 = seq(0, 1.25, length.out = 1043)
#' random_ratio = 1
#'
#'
#' plot(out)
#'
#' class(p)
#' typeof(p)
#' p
#' power_gsnb(ratio_H1 = 0.7,
#'            timing = c(0.5, 1), power_gs = 0.8, sig_level = 0.025)
#'
#' x <- 1:10
#' y <- x + rnorm(length(x))
#' lreg <- lm(y~x)
#' plot(lreg)
#'
#' timing <- c(0.5, 1)
#' sig_level <- 0.025
#' max_info <- 61.93
#' esf_out <- esf_obrien(t = timing, sig_level = sig_level)
#' alloc_error <- c(esf_out[1], diff(esf_out))
#' ## Calculate the critical values
#' critical <- get_critical_values(timing = timing, alloc_error = alloc_error)
#' get_rejectprob_gsnb(rate_ratio = 1, ratio_H0 = 1, critical = critical, max_info, timing)
#'
#'
#'
#'
#' library(gsDesign)
#' gsDesign(k=2, test.type=1, delta=-log(0.7), sfu = "Pocock")
#' out$critical
#'
#' f <- function() {
#'   a <- 1
#'  body_a <- quote({a <- 2})
#'  eval(body_a)
#'   a
#' }
#' f()
#
#
#
# Calculate study period for given recruitment times
# library(gscounts)
# expose <- seq(0, 1.25, length.out = 1042)
# out <- design_nb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power = 0.8,
#                    ratio_H0 = 1, sig_level = 0.025,
#                    t_recruit1 = expose, t_recruit2 = expose, random_ratio = 1)
# out
# # Calculate sample size for a fixed exposure time
# out <- design_nb(rate1 = 4.2, rate2 = 8.4, shape = 3, power = 0.8,
#                    ratio_H0 = 1, sig_level = 0.025,
#                    followup_max = 0.5, random_ratio = 1)
# out$n1
#  (qnorm(1-0.025) + qnorm(0.8))^2 / log(0.5)^2 * ((4.2+8.4) / (0.5*4.2*8.4) +2 * 3)

#
#
#
# # testing the critical values
# library("ldbounds")
# library("gsDesign")
#
# esf_obrien(t = c(0.5, 1), sig_level = 0.025)
# system.time({
#   for(i in 1:1000) {
#     get_critical_values(timing = c(0.5, 1), alloc_error = c(0.001525323, 0.025-0.001525323))
#   }
# })
# system.time({
#   for(i in 1:1000) {
#     bounds(t = c(0.5, 1), iuse=1,alpha=0.025)
#   }
# })
# system.time({
#   for(i in 1:1000) {
#     gsDesign(k = 2, test.type=2, sfu=sfLDOF, alpha=.025)$upper$bound
#   }
# })
#
# gsDesign(k = 2, test.type=2, sfu=sfLDOF, alpha=.025)
# bounds(t = c(0.5, 1), iuse=1,alpha=0.025)$upper.bounds
# gsDesign(k = 2, test.type=1, sfu=sfLDOF, alpha=.025)$upper$bound
#
#
#
# esf_obrien(t=c(0.5, 1), sig_level = 0.05)
# esf_obrien2(t=c(0.5, 1), sig_level = 0.025)
#
# esf_obrien2 <- function(t, sig_level, ...) {
#   if (any(t<0)) stop("t must be non-negative")
#   f <- 1- 1 * pnorm(qnorm(1-sig_level)/sqrt(t))
#   pmin(f, sig_level)
# }
#
# upper.bounds
# bounds(t = 1, iuse=1,alpha=0.025-0.001525323)$upper.bounds
# bounds(t = c(1), iuse=1,alpha=0.025-0.001525323)$upper.bounds
