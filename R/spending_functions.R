#' pocock
#' @description Error spending function mimicking Pococks critical values
#' @param t numeric; Non-negative information ratio
#' @param sig_level numeric; significance level
#' @param ... optional arguments
#' @return numeric
#' @examples 
#' # Pocock-type error spending function
#' pocock(t = c(0.5, 1), sig_level = 0.025)
#' @export
pocock <- function(t, sig_level, ...) {
  if (any(t < 0)) 
    stop("t must be non-negative")
  f <- sig_level * log(1 + (exp(1) - 1) * t)
  pmin(f, sig_level)
}


#' obrien
#' @description Error spending function mimicking O'Brien & Fleming critical values
#' @param t numeric; Non-negative information ratio
#' @param sig_level numeric; significance level
#' @param ... optional arguments
#' @return numeric
#' @import stats
#' @examples 
#' # O'Brien-Fleming-type error spending function
#' obrien(t = c(0.5, 1), sig_level = 0.025)
#' @export
obrien <- function(t, sig_level, ...) {
  if (any(t < 0)) 
    stop("t must be non-negative")
  f <- 2 - 2 * pnorm(qnorm(1-sig_level/2) / sqrt(t))
  pmin(f, sig_level)
}
