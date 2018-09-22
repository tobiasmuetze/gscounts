# context("Test various recruitment times")
# 
# 
# design1 <- design_gsnb(rate1 = 1, rate2 = 1.5, dispersion = 3,
#                    power = 0.8, timing = c(0.5, 0.75, 1), esf = pocock,
#                    ratio_H0 = 1, sig_level = 0.025, 
#                    study_period = 2,
#                    accrual_period = 1, random_ratio = 1)
# design1
# 
# design2 <- design_gsnb(rate1 = 1, rate2 = 1.5, dispersion = 3,
#                        power = 0.8, timing = c(0.5, 0.75, 1), 
#                        esf = pocock,
#                        ratio_H0 = 1, sig_level = 0.025, 
#                        t_recruit1 = seq(0, 1, length.out = 399), 
#                        t_recruit2 = seq(0, 1, length.out = 399), 
#                        random_ratio = 1)
# design2
# 
# 
# design3 <- design_gsnb(rate1 = 1, rate2 = 1.5, dispersion = 3,
#                        power = 0.8, timing = c(0.5, 0.75, 1), 
#                        esf = pocock,
#                        ratio_H0 = 1, sig_level = 0.025, 
#                        study_period = 2,
#                        accrual_period = 1,
#                        accrual_speed = 2,
#                        random_ratio = 1)
# design3
# 
# 
# design1$calendar
# design2$calendar
# design3$calendar
