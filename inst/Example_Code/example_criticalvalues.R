rm(list = ls())


get_critical_values(timing = x$timing, alloc_error = x$efficacy$spend)
x$efficacy$critical

cpp_calc_critical(r = 20, lower = -0.2151421, upper = 3.356869, error_spend = 0.05661261, information = c(4.571987, 9.143974), theta = 1, side = "lower")

.Call('gscounts_cpp_calc_critical', PACKAGE = 'gscounts', 
      r = 20, lower = -0.2151421, upper = 3.356869, error_spend = 0.05661261, 
      information = c(4.571987, 9.143974), theta = 1, side = "lower")


gscounts::get_critical_values(timing = c(0.5, 0.75, 1), alloc_error = c(0.01, 0.02, 0.02))

.Call('gscounts_cpp_calc_critical', PACKAGE = 'gscounts', 
      r = 20, lower = -2.326348, upper = Inf, error_spend = 0.02, 
      information = c(1, 1.5), theta = 0, side = "lower")


# Critical values without futility
get_critical_values(timing = c(0.2, 0.5, 1), alloc_error = c(1.172645e-05, 0.00556287, 0.04442540))
get_critical_values2(timing = c(0.2, 0.5, 1), error_futility = NULL, error_efficacy = c(1.172645e-05, 0.00556287, 0.04442540), futility = "none")


# Critical values with non-binding futility
get_critical_values2(timing = c(0.2, 0.5, 1), error_futility = c(1.172645e-05, 0.00556287, 0.04442540), 
                     error_efficacy = c(1.172645e-05, 0.00556287, 0.04442540), 
                     futility = "nonbinding", log_effect = -1)

get_critical_values(timing = c(0.2, 0.4, 0.8), alloc_error = c(0.01, 0.02, 0.02))

library(gsDesign)
sfLDOF(alpha = 0.05, t = c(0.2, 0.5, 1))$spend
diff(sfLDOF(alpha = 0.05, t = c(0.2, 0.5, 1))$spend)

x <- gsDesign(k=3, timing = c(0.2, 0.5, 1), test.type=4, sfu = sfLDOF, sfl = sfLDOF, delta = 1, alpha = 0.05, beta = 0.05)
plot(x)
x$lower$bound
x$upper$bound

x$lower$spend
x$upper$spend

qnorm(x$upper$spend[1], mean = sqrt(0.2) )

system.time({
  for(i in 1:1000)
    get_critical_values(timing = c(0.2, 0.4, 0.8), alloc_error = c(0.01, 0.02, 0.02))
})

system.time({
  for(i in 1:1000)
    get_critical_values2(timing = c(0.2, 0.4, 0.8), error_futility = NULL, error_efficacy = c(0.01, 0.02, 0.02), futility = "none")
})

