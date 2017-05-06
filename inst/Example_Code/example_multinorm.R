library(mvtnorm)

a1 <- -Inf
b1 <- 1.5
a2 <- -Inf
b2 <- 0.8
a3 <- -Inf
b3 <- 3
a4 <- -2
b4 <- 5
a5 <- -1
b5 <- Inf

pmvnorm(lower = c(a1, a2, a3, a4, a5), upper = c(b1, b2, b3, b4, b5), 
        sigma = get_covar(timing = c(0.2, 0.4, 0.6, 0.8, 1)))
pmultinorm(r = 18, lower = c(a1, a2, a3, a4, a5), upper = c(b1, b2, b3, b4, b5), 
           information = c(0.2, 0.4, 0.6, 0.8, 1), theta = 0)

.Call('gscounts_cpp_pmultinorm', PACKAGE = 'gscounts', 
      r = 20, lower = c(a1, a2), upper = c(b1, b2), information = c(0.4, 0.6), theta = 1)
pmvnorm(lower = c(a1, a2), upper = c(b1, b2), 
        sigma = get_covar(timing = c(0.4, 0.6)), mean = sqrt(c(0.4, 0.6)))


pmultinorm(r = 18, lower = c(a1, a2), upper = c(b1, b2), information = c(0.4, 0.6), theta = 0)
pmultinorm(r = 18, lower = c(b1, b2), upper = c(Inf, Inf), information = c(0.4, 0.6), theta = 0)



system.time({
  for(i in 1:1000) {
    pmultinorm(r = 10, lower = c(a1, a2, a3), upper = c(b1, b2, b3), information = c(0.25, 0.5, 1), theta = 0)
  }
})
system.time({
  for(i in 1:1000) {
    pmvnorm(lower = c(a1, a2, a3), upper = c(b1, b2, b3), sigma = get_covar(timing = c(0.25, 0.5, 1)))
  }
})