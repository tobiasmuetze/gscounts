#' @name check_gsnb
#' @title Check object of class gsnb
#' @description Performs validity checks of object of class gsnb
#' @param x object of class gsnb
#' @keywords internal
check_gsnb <- function(x) {
  
  if (class(x) != "gsnb")
    stop("x ist not an object of class gsnb")
  
  # rate1
  if (is.null(x$rate1))
    stop("rate1 is not defined")
  if ((x$rate1 <= 0))
    stop("rate1 must be positive")
  # rate2
  if (is.null(x$rate2))
    stop("rate2 is not defined")
  if ((x$rate2 <= 0))
    stop("rate2 must be positive")
  # dispersion
  if (is.null(x$dispersion))
    stop("dispersion is not defined")
  if ((x$dispersion <= 0))
    stop("dispersion must be positive")
  # rate_H0
  if (is.null(x$ratio_H0))
    stop("rate_H0 is not defined")
  if ((x$ratio_H0 <= x$rate1 / x$rate2))
    stop("rate1/rate2 is not located in the alternative")
  # significance level
  if (is.null(x$sig_level))
    stop("Power is not defined")
  if ((x$sig_level <= 0) || (x$sig_level >= 1))
    stop("sig_level must be between 0 and 1")
  # power
  if (is.null(x$power))
    stop("Power is not defined")
  if ((x$power <= x$sig_level) || (x$power >= 1))
    stop("power must be greater than significance level and smaller than 1")
  # efficacy
  if (is.null(x$efficacy))
    stop("efficacy parameters are not specified")
  if (!is.function(x$efficacy$esf))
    stop("esf is not an error spending function")
  if (abs(sum(x$efficacy$spend) - x$sig_level) > 1e-6)
    stop("efficacy spending does not sum up to significance level")
  # futility
  if (!is.null(x$futility)) {
    if (!(x$futility$type %in% c("binding", "nonbinding")))
      stop("futility must be either binding or non-binding")
    if (!is.function(x$futility$esf))
      stop("esf_futility is not an error spending function")
    if (abs(sum(x$futility$spend) - 1 + x$power) > 1e-6)
      stop("futility spending does not sum up to 1-power")
  }
}

