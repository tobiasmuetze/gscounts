#' @name get_critical_values
#' @title Calculate critical values
#' @description Calculate either critical value for analysis at stage k based on
#' given previous values or calculate set of critical values for multiple analysis.
#' @details If \code{alloc_error} has length 1 and \code{prev_critical} is \code{NULL},
#' the corresponding normal quantile, i.e. the critical value for the first interim,
#' is returned.
#' If length of \code{alloc_error} and \code{timing} is identical and \code{prev_critical}
#' is \code{NULL}, the vector of critical values for the analyses at times \code{timing}
#' is returned.
#' If \code{alloc_error} has length 1 and \code{timing} has one entry more than
#' \code{prev_critical}, the critical value for the analysis at time \code{timing[k]}
#' is returned.
#' @param timing numeric vector; 0 < \code{timing[1]} < ... < \code{timing[k]} with
#' \code{k} the index and \code{timing[k]} the information
#' ratio of the current analysis
#' @param alloc_error Numeric vector with allocated type I error.
#' @param prev_critical Critical values of previous analyses
#' @return numeric; critical values
#' @export
#' @import mvtnorm
get_critical_values <- function(timing, alloc_error, prev_critical = NULL) {

  k <- length(timing)
  # Calculate the canonical form covariance matrix
  covar <- get_covar(timing)
  # Calculation of rejection probability for stage i
  calc_rejectprob <- quote({
    lower <- c(critical, -Inf)
    upper <- c(rep(Inf, times = i-1), x)
    pmvnorm(lower = lower, upper = upper, mean = rep(0, times = i), sigma = covar[1:i, 1:i])
  })

  if (is.null(prev_critical) && (length(timing) == length(alloc_error)) ) {
    # Calculate critical values for each entry of alloc_error
    critical <- qnorm(alloc_error[1])
    for(i in 2:k) {
      critical[i] <- uniroot(f = function(x) eval(calc_rejectprob) - alloc_error[i],
                             interval = c(-10, -1e-6), #extendInt = "downX",
                             tol = .Machine$double.eps^0.5)$root
    }
    return(critical)
  } else if (!is.null(prev_critical) && ((k-1) == length(prev_critical)) && (length(alloc_error) == 1)) {
    if (alloc_error == 0) {
      return(-Inf)
    }
    # Find quantile such that the allocated error is spend
    i <- k
    critical <- prev_critical
    c_k <- uniroot(f = function(x) eval(calc_rejectprob) - alloc_error,
                   interval = c(-10, -1e-6), #extendInt = "downX",
                   tol = .Machine$double.eps^0.5)$root
    return(c_k)
  } else if (is.null(prev_critical) && (length(alloc_error) == 1)) {
    return(qnorm(alloc_error))
  } else {
    stop("Wrong arguments")
  }

}
