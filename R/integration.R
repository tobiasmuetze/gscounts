#' rm(list = ls())
#' library(mvtnorm)
#'
#' #' get_covar
#' #' @description Get canonical form covariance matrix
#' #' @param timing numeric vector; 0 < \code{timing[1]} < ... < \code{timing[k]}
#' #' @keywords internal
#' get_covar <- function(timing) {
#'   k <- length(timing)
#'   # Calculate matrix with (i,j) = I_i / I_j
#'   covar <- sweep(x = matrix(timing, nrow = k, ncol = k), MARGIN = 2, STATS = timing, FUN = '/')
#'   # Make matrix symmetric
#'   covar[lower.tri(covar)] <- t(covar)[lower.tri(covar)]
#'   covar <- sqrt(covar)
#'   covar
#' }
#'
#'
#' get_gridpoints <- function(r, lower, upper, add_points = NULL) {
#'
#'   # m <- 12 * r - 3
#'   # gridpoints <- weight <- numeric(m)
#'   #
#'
#'   x <- -3 - 4 * log(r / (1:(r-1)))
#'   x <- c(x, -3 + 3 * (0:(4*r)) / (2*r))
#'   x <- c(x, 3 + 4*log(r/((r-1):1)))
#'   x <- sort(c(x, add_points))
#'
#'   if (any(x < lower)) {
#'     x <- x[x > lower]
#'     x <- x[x > lower]
#'     x <- c(lower, x)
#'   }
#'   if (any(x > upper)) {
#'     x <- x[x > lower]
#'     x <- x[x > lower]
#'     x <- c(x, upper)
#'   }
#'
#'   m <- 2*length(x) - 1
#'   gridpoints <- weight <- numeric(m)
#'
#'   odd_indexes <- seq(1, m, by = 2)
#'   even_indexes <- seq(2, m, by = 2)
#'
#'   gridpoints[odd_indexes] <- x
#'   gridpoints[even_indexes] <- (head(x, -1) + tail(x, -1)) / 2
#'
#'   # gridpoints <- gridpoints[gridpoints < lower]
#'   # gridpoints <- c(gridpoints, lower)
#'   # m <- length(gridpoints)
#'   # weight <- numeric(m)
#'
#'   odd_indexes <- seq(3, m-2, by = 2)
#'   even_indexes <- seq(2, m-1, by = 2)
#'   weight[1] <- (gridpoints[3] - gridpoints[1]) / 6
#'   weight[m] <- (gridpoints[m] - gridpoints[m-2]) / 6
#'   weight[odd_indexes] <- (gridpoints[odd_indexes+2] - gridpoints[odd_indexes-2]) / 6
#'   weight[even_indexes] <- 4 * (gridpoints[even_indexes+1] - gridpoints[even_indexes-1]) / 6
#'   list(gridpoints = gridpoints, weight = weight)
#' }
#'
#' get_gridpoints(r = 2, lower = -2, upper = 2)
#'
#'
#' f_k <- function(zk, zk1, theta, Ik, Ik1) {
#'   Delta <- Ik - Ik1
#'
#'   x <- (zk * sqrt(Ik) - zk1 * sqrt(Ik1) - theta * Delta) / sqrt(Delta)
#'   sqrt(Ik / Delta) * dnorm(x = x, mean = 0, sd = 1)
#' }
#'
#' e_k1 <- function(zk1, Ik1, theta, bk, Ik) {
#'   Delta <- Ik - Ik1
#'   x <- (zk1 * sqrt(Ik1) - bk * sqrt(Ik) + theta * Delta) / sqrt(Delta)
#'   pnorm(x, mean = 0, sd = 1)
#' }
#'
#' update_h <- function(h, gridptsIk1, lower, upper, r, theta, Ik, Ik1) {
#'
#'   gridpoint_out <- get_gridpoints(r = r, lower = lower, upper = upper)
#'   gridpnt_new <- gridpoint_out$gridpoints
#'   wght <- gridpoint_out$weight
#'
#'   h_new <- numeric(length(gridpnt_new))
#'   for(i in seq_along(gridpnt_new)) {
#'     h_new[i] <- sum(f_k(zk = gridpnt_new[i], zk1 = gridptsIk1, theta = theta, Ik = Ik, Ik1 = Ik1) * wght[i] * h)
#'   }
#'   list(h = h_new, weight = wght, gridpnts = gridpnt_new)
#' }
#'
#'
#'
#' a1 <- -2.09
#' b1 <- 2.08
#' a2 <- -1.789
#' b2 <- 1.684
#' a3 <- 1.5
#' b3 <- Inf
#'
#' update_h_out <- update_h(h = 1, gridptsIk1 = 0, lower = a1, upper = b1, r = 18, theta = 0, Ik = 30, Ik1 = 0)
#' update_h_out2 <- update_h(h = update_h_out$h, gridptsIk1 = update_h_out$gridpnts, lower = a2, upper = b2, r = 18, theta = 0, Ik = 60, Ik1 = 30)
#' sum(update_h_out2$h * e_k1(zk1 = update_h_out2$gridpnts, Ik1 = 60, theta = 0, bk = a3, Ik = 90))
#'
#'
#' update_h_out <- update_h(h = 1, gridptsIk1 = 0, lower = a1, upper = b1, r = 18, theta = 0, Ik = 1/3, Ik1 = 0)
#' update_h_out2 <- update_h(h = update_h_out$h, gridptsIk1 = update_h_out$gridpnts, lower = a2, upper = b2, r = 18, theta = 0, Ik = 2/3, Ik1 = 1/3)
#' update_h_out3 <- update_h(h = update_h_out2$h, gridptsIk1 = update_h_out2$gridpnts, lower = a3, upper = b3, r = 18, theta = 0, Ik = 1, Ik1 = 2/3)
#' sum(update_h_out3$h)
#'
#' pmvnorm(lower = c(a1, a2, a3), upper = c(b1, b2, b3), sigma = get_covar(timing = c(1/3, 2/3, 1)), algorithm = "Miwa")[1]
