
# sink("/home/tobiasmuetze/Dropbox/NVS Internship/R-Package/gscounts/inst/Example_Code/output_design_gsnb.txt")


# Calculate the sample sizes for a given accrual period and study period
out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
                   timing = c(0.5, 1), esf = esf_obrien,
                   ratio_H0 = 1, sig_level = 0.025,
                   study_period = 3.5, accrual_period = 1.25, random_ratio = 1)
out

# Calculate study period for given recruitment times
expose <- seq(0, 1.25, length.out = 1042)
out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
                   timing = c(0.5, 1), esf = esf_obrien,
                   ratio_H0 = 1, sig_level = 0.025,
                   t_recruit1 = expose, t_recruit2 = expose, random_ratio = 1)
out

# Calculate sample size for a fixed exposure time
out <- design_gsnb(rate1 = 0.0875, rate2 = 0.125, shape = 5, power_gs = 0.8,
                   timing = c(0.5, 1), esf = esf_obrien,
                   ratio_H0 = 1, sig_level = 0.025,
                   followup_max = 0.5, random_ratio = 1)
out
#sink()
