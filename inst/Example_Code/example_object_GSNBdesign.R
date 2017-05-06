rm(list = "x")
x <- list(rate1 = 0.0875, rate2 = 0.125, 
          dispersion = 5, ratio_H0 = 1,
          power = 0.8, sig_level = 0.05,
          timing = c(0.2, 0.5, 1), esf = esf_obrien,
          esf_futility = esf_obrien, 
          followup_max = 0.5, random_ratio = 1)
x$efficacy$esf <- esf_obrien
x$efficacy$spend <- c(1.172645e-05, diff(esf_obrien(t = x$timing, sig_level = 0.05)))
x$futility$spend <- c(0.004161719, diff(esf_obrien(t = x$timing, sig_level = 0.2)))
x$futility$type <- "binding"
x <- add_maxinfo(x)
x$efficacy$critical
x$futility$critical

library(gsDesign)
# Two-sided test with binding futility boundaries
y <- gsDesign(k=3, timing = c(0.2, 0.5, 1), test.type=3, sfu = sfLDOF, sfl = sfLDOF, delta = -log(0.7), alpha = 0.05, beta = 0.2)
y$lower$bound
y$upper$bound
y
plot(y)


# Two-sided test with non-binding futility boundaries
y <- gsDesign(k=3, timing = c(0.2, 0.5, 1), test.type=4, sfu = sfLDOF, sfl = sfLDOF, delta = -log(0.7), alpha = 0.05, beta = 0.2)
y$lower$bound
y$upper$bound
y
plot(y)



# One-sided test without futility boundaries
y <- gsDesign(k=3, timing = c(0.2, 0.5, 1), test.type=1, sfu = sfLDOF, sfl = sfLDOF, delta = -log(0.7), alpha = 0.05, beta = 0.2)
y$lower$spend
y$upper$bound
y
plot(y)
power_gsnb(ratio_H1 = 0.7, power_gs = 0.8, timing = x$timing, esf = esf_obrien, ratio_H0 = 1, sig_level = 0.05)
get_critical_values(timing = c(0.2, 0.5, 1), alloc_error = c(8.857544e-05, 5.486021e-03, 4.442540e-02))

qnorm(0.004161719, mean = 0.3566749 * sqrt(0.2*52))


get_covar <- function(timing) {
  k <- length(timing)
  # Calculate matrix with (i,j) = I_i / I_j
  covar <- sweep(x = matrix(timing, nrow = k, ncol = k), MARGIN = 2, STATS = timing, FUN = '/')
  # Make matrix symmetric
  covar[lower.tri(covar)] <- t(covar)[lower.tri(covar)]
  covar <- sqrt(covar)
  covar
}
y <- gsDesign(k=3, timing = c(0.2, 0.5, 1), test.type=1, sfu = sfLDOF, alpha = 0.05, beta = 0.2, delta = 0.5)
plot(y)
y
1-pmvnorm(lower = rep(-Inf, times = 3), upper = y$upper$bound, sigma = get_covar(timing = c(0.2, 0.5, 1)))
1-pmvnorm(lower = rep(-Inf, times = 3), upper = y$upper$bound, sigma = get_covar(timing = c(0.2, 0.5, 1)), mean = 0.5*sqrt(c(0.2, 0.5, 1) * 48))

sfup <- c(.033333, .063367, .1)
sflp <- c(.25, .5, .75)
timing <- c(.1, .4, .7, 1)
x <- gsDesign(k=4, timing=timing, sfu=sfPoints, sfupar=sfup, sfl=sfPoints,
              sflpar=sflp,n.fix=1000) 
x
plot(x)
1-pmvnorm(lower = rep(-Inf, times = length(timing)), upper = x$upper$spend, sigma = get_covar(timing = timing))
