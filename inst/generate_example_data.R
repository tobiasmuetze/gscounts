rm(list = ls())

# Generate events for group 1 and 2
poisson_rates1 <- rgamma(n = out$n1, shape = 1/dispersion, rate = 1/(dispersion * rate1))
poisson_rates2 <- rgamma(n = out$n2, shape = 1/dispersion, rate = 1/(dispersion * rate2))
hospitalizations <- data.frame()
## Generate events for group 1 / experimental group
for(i in seq_along(poisson_rates1)) {
  counts <- rpois(n = 1, lambda = (study_period - out$t_recruit1[i]) * poisson_rates1[i])
  if (counts == 0) {
    hospitalizations <- rbind(hospitalizations, 
                              data.frame(treatment = "experiment",
                                         pat = i,
                                         t_recruit = out$t_recruit1[i],
                                         eventtime = NA))
  } else {
    hospitalizations <- rbind(hospitalizations, 
                              data.frame(treatment = "experiment",
                                         pat = i,
                                         t_recruit = out$t_recruit1[i],
                                         eventtime = runif(n = counts, min = out$t_recruit1[i], 
                                                           max = study_period)))
  }
}

## Generate events for group 2 / control group
for(i in seq_along(poisson_rates2)) {
  counts <- rpois(n = 1, lambda = (study_period - out$t_recruit2[i]) * poisson_rates2[i])
  if (counts == 0) {
    hospitalizations <- rbind(hospitalizations, 
                              data.frame(treatment = "control",
                                         pat = i,
                                         t_recruit = out$t_recruit2[i],
                                         eventtime = NA))
  } else {
    hospitalizations <- rbind(hospitalizations, 
                              data.frame(treatment = "control",
                                         pat = i,
                                         t_recruit = out$t_recruit2[i],
                                         eventtime = runif(n = counts, min = out$t_recruit2[i], 
                                                           max = study_period)))
  }
}



monitor_info <- sapply(X = tail(sort(unique(hospitalizations$eventtime)), -5), 
                       FUN = function(t) {
                         x <- filter(hospitalizations, t_recruit <= t, (is.na(eventtime) | eventtime <= t))
                         x <- x %>% 
                           group_by(treatment, pat, t_recruit) %>% 
                           summarise(count = sum(eventtime <= t, na.rm = TRUE)) 
                         x$t_expose <- t - x$t_recruit
                         glm_out <- glm.nb(data = x, formula = count ~ treatment + offset(log(t_expose)) - 1)
                         rates <- exp(glm_out$coefficients)
                         followup1 <- filter(x, treatment == 'experiment')$t_expose
                         followup2 <- filter(x, treatment == 'control')$t_expose
                         
                         get_info_gsnb(rate1 = rates["treatmentexperiment"],
                                       rate2 = rates["treatmentcontrol"],
                                       dispersion = 1/glm_out$theta, 
                                       followup1 = followup1, 
                                       followup2 = followup2)
                       })

