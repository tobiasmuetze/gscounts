rm(list = ls())

# Calculate the sample sizes for a given accrual period and study period
study_period <- 4
rate1 <- 0.0875
rate2 <- 0.125
dispersion <- 5
out <- design_gsnb(rate1 = rate1, rate2 = rate2, dispersion = dispersion, 
                   power = 0.8,
                   timing = c(0.4, 0.7, 1), esf = obrien,
                   ratio_H0 = 1, sig_level = 0.025,
                   study_period = study_period, accrual_period = 1.25, random_ratio = 1)
out


rm(list = ls())

out4 <- design_gsnb(rate1 = 4.2, rate2 = 8.4, dispersion = 3, ratio_H0 = 1, 
                    power = 0.8, sig_level = 0.025, timing = c(0.5, 1), 
                    esf = pocock, followup_max = 0.5, random_ratio = 2,
                    futility = "nonbinding", esf_futility = obrien)

out4


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

#devtools::use_data(hospitalizations, hospitalizations)

library(gscounts)
library(dplyr)
library(MASS)
x <- hospitalizations %>% 
  group_by(treatment, pat) %>% 
  summarise(count = sum(na.rm = TRUE))
x


x <- hospitalizations %>% 
  group_by(treatment, pat) %>% 
  summarise(count = sum(!is.na(eventtime))) 
x

x <- hospitalizations %>% 
  group_by(treatment, pat, t_recruit) %>% 
  summarise(count = sum(eventtime < 0.551, na.rm = TRUE)) 
x

head(hospitalizations)

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

plot(tail(sort(unique(hospitalizations$eventtime)), -5), monitor_info)
tail(sort(unique(hospitalizations$eventtime)), -5)[123]
which(monitor_info > 62.66371*0.4)

t <- 1.297192
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

Rcpp::sourceCpp(file = "/home/tmuetze/Dropbox/NVS Internship/Statistical Methods in Medical Research/R/Cpp/estimate_nb.cpp")
estimate_nb(counts1 = filter(x, treatment == 'experiment')$count, 
            counts2 = filter(x, treatment == 'control')$count, 
            t1 = filter(x, treatment == 'experiment')$t_expose, 
            t2 = filter(x, treatment == 'control')$t_expose)
glm_out <- glm.nb(data = x, formula = count ~ treatment + offset(log(t_expose)) )
summary(glm_out)
1/sqrt(25.36104)
# %>% 
#   count(count)
