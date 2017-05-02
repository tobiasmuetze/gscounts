library(gscounts)

########################################################
# Study period for given recruitment time
########################################################
t_recruit1 <- seq(0, 1.25, length.out = 1200)
t_recruit2 <- seq(0, 1.25, length.out = 800)
out <- design_nb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power = 0.8,
                 ratio_H0 = 1, sig_level = 0.025, random_ratio = 1.5,
                 t_recruit1 = t_recruit1, t_recruit2 = t_recruit2)
out

study_period <- 4.020687
exposure1 <- study_period - t_recruit1
exposure2 <- study_period - t_recruit2

pvalue <- numeric(1000)
for(i in seq_along(pvalue)) {
  x1 <- rnbinom(n = length(t_recruit1), mu = exposure1 * 0.0875, size = 1 / 5)
  x2 <- rnbinom(n = length(t_recruit2), mu = exposure2 * 0.125, size = 1 / 5)
  group <- c(rep("A", times = length(t_recruit1)), rep("B", times = length(t_recruit2)))
  nb_fit <- MASS::glm.nb(c(x1, x2) ~ group + offset(c(exposure1, exposure2)))
  pvalue[i] <- 1-pnorm(summary(nb_fit)$coeff[2, 3]) # pvalue for test of H0:rate1 >= rate2
}
mean(pvalue <= 0.025)
########################################################
# END Study period for given recruitment time
########################################################




########################################################
# Sample size for a fixed exposure time
########################################################

# Calculate sample size for given accrual period and study duration assuming uniformal accrual
out <- design_nb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power = 0.8,
                 ratio_H0 = 1, sig_level = 0.025,
                 study_period = 4, accrual_period = 1, random_ratio = 2)
out

pvalue <- numeric(1000)
for(i in seq_along(pvalue)) {
  x1 <- rnbinom(n = out$n1, mu = 0.5 * 4.2, size = 1 / 3)
  x2 <- rnbinom(n = out$n2, mu = 0.5 * 8.4, size = 1 / 3)
  group <- c(rep("A", times = out$n1), rep("B", times = out$n2))
  nb_fit <- MASS::glm.nb(c(x1, x2) ~ group + offset(c(rep(0.5, times = out$n1), rep(0.5, times = out$n2))))
  pvalue[i] <- 1-pnorm(summary(nb_fit)$coeff[2, 3])
}
mean(pvalue <= 0.025)

########################################################
# END Sample size for a fixed exposure time
########################################################


#####################################################################################
# Sample size for given accrual period and study duration assuming uniformal accrual
#####################################################################################
out <- design_nb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power = 0.8,
                 ratio_H0 = 1, sig_level = 0.025,
                 study_period = 4, accrual_period = 1, random_ratio = 2)
out

exposure1 <- out$study_period - out$t_recruit1
exposure2 <- out$study_period - out$t_recruit2

pvalue <- numeric(2000)
for(i in seq_along(pvalue)) {
  x1 <- rnbinom(n = length(exposure1), mu = exposure1 * 0.0875, size = 1 / 5)
  x2 <- rnbinom(n = length(exposure2), mu = exposure2 * 0.125, size = 1 / 5)
  group <- c(rep("A", times = length(exposure1)), rep("B", times = length(exposure2)))
  nb_fit <- MASS::glm.nb(c(x1, x2) ~ group + offset(c(exposure1, exposure2)))
  pvalue[i] <- 1-pnorm(summary(nb_fit)$coeff[2, 3]) # pvalue for test of H0:rate1 >= rate2
}
mean(pvalue <= 0.025)
#####################################################################################
# END Sample size for given accrual period and study duration assuming uniformal accrual
###################################################################################