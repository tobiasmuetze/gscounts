#' power_gsnb
#' @description Calculate power of group sequential design with
#' negative binomial endpoints
#' @param ratio_H0 numeric; postive number denoting the rate ratio rate_1/rate_2
#' under the null hypothesis, i.e. the non-inferiority or superiority margin
#' @param ratio_H1 numeric; postive number denoting the rate ratio rate_1/rate_2
#' under the alternative hypothesis
#' @param max_info The maximum information of the study.
#' @param timing numeric vector; 0 < timing[1] < ... < timing[K] = 1 with K the
#' number of analyses, i.e. (K-1) interim analyses and final analysis
#' @param esf An error spending function
#' @param power_gs numeric; the target power of the group sequential design
#' @param sig_level numeric; the total significance level
#' @param ... additional arguments for the error spending function esf
#' @return A list with class "GSNBpower" containing the following components:
#' \item{ratio_H0}{as input}
#' \item{ratio_H1}{as input}
#' \item{timing}{as input}
#' \item{sig_level}{as input}
#' \item{max_info}{as input if not NULL; otherweise calculated maximum information}
#' \item{power_gs}{as input if not NULL; otherweise calculated power}
#' \item{power_fix}{power of fixed design}
#' \item{critical}{critical values for each look}
#' @import mvtnorm
#' @import stats
#' @export
power_gsnb <- function(ratio_H1, max_info = NULL, power_gs = NULL, timing, esf = esf_obrien, ratio_H0 = 1, sig_level, ...) {
  
  # Parameter checks
  stopifnot((ratio_H1 > 0) && (ratio_H1 < ratio_H0))
  stopifnot(all(timing > 0) && all(timing <= 1) && (tail(timing, n = 1) == 1) && all(diff(timing) > 0))
  stopifnot(is.null(max_info) || (max_info > 0))
  stopifnot((sig_level > 0) && (sig_level < 1))
  stopifnot(is.null(power_gs) || ((power_gs > sig_level) && (power_gs < 1)))
  
  arguments <- c(as.list(environment()), list(...))
  
  # Number of analyses and effect size
  K <- length(timing)
  effect_size <- log(ratio_H1) - log(ratio_H0)
  
  # Calculate the canonical form covariance matrix
  covar <- get_covar(timing)
  
  # Critical values
  ## Vector with type I errors for each stage
  esf_out <- esf(t = timing, sig_level = sig_level, ...)
  alloc_error <- c(esf_out[1], diff(esf_out))
  ## Calculate the critical values
  critical <- get_critical_values(timing = timing, alloc_error = alloc_error)
  
  # Power in group sequential design
  power_body <- quote({
    info_vec <- max_info * timing
    mean_vec <- effect_size * sqrt(info_vec)
    1 - pmvnorm(lower = critical, upper = Inf, mean = mean_vec, sigma = covar)[1]
  })
  
  # Calculate the missing value
  if (is.null(power_gs)) {
    power_gs <- eval(power_body)
    power_fix <- pnorm(qnorm(sig_level), mean = effect_size * sqrt(max_info))
  } else if (is.null(max_info)) {
    max_info <- uniroot(f = function(max_info){eval(power_body) - power_gs},
                        interval = c(0.1, 10), extendInt = "upX")$root
    power_gs <- eval(power_body)
    power_fix <- pnorm(qnorm(sig_level), mean = effect_size * sqrt(max_info))
  }


  out <- c(arguments[c("ratio_H1", "timing", "ratio_H0", "sig_level", "esf")],
           power_gs = power_gs, power_fix = power_fix, max_info = max_info, list(critical = critical))
  class(out) <- "GSNBpower"
  out
}
