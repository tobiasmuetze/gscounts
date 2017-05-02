#' print.GSNBpower
#' @description print method for instance of class GSNBpower
#' @param x an object of class GSNBpower
#' @param ... optional arguments to print or plot methods
#' @export
print.GSNBpower <- function(x, ...) {
  cat("\t Power for group sequential design with negative binomial outcomes\n\n")
  cat("Rate ratio under null hypothesis: ", format(x$ratio_H0, digits = 4), "\n")
  cat("Rate ratio under alternative: ", format(x$ratio_H1, digits = 4), "\n")
  cat("Significance level: ", format(x$sig_level, digits = 4), "\n")
  cat("Analyses: ", length(x$timing), "\n")
  cat("Information times: ", paste(format(x$timing, digits = 4), collapse = ", "), "\n")
  cat("Maximum information: ", paste(format(x$max_info, digits = 4), collapse = ", "), "\n")
  cat("Power group sequential design: ", format(x$power_gs, digits = 4), "\n")
  cat("Power fixed design: ", format(x$power_fix, digits = 4), "\n")
}


#' print.GSNBdesign
#' @description print method for instance of class GSNBpower
#' @param x an object of class GSNBpower
#' @param ... optional arguments to print or plot methods
#' @export
print.GSNBdesign <- function(x, ...) {
  cat("\t Design of group sequential trial with negative binomial outcomes\n\n")
  cat("Rate ratio under null hypothesis: ", format(x$ratio_H0, digits = 4), "\n")
  cat("Rate ratio under alternative: ", format(x$rate1 / x$rate2, digits = 4), "\n")
  cat("Rate group 1: ", x$rate1, "\n")
  cat("Rate group 2: ", x$rate2, "\n")
  cat("Shape parameter: ", x$shape, "\n")
  cat("Number of looks: ", length(x$timing), "\n")
  if (!is.null(x$accrual_period)) cat("Accrual period: ", x$accrual_period, "\n")
  if (!is.null(x$study_period)) cat("Study duration: ", x$study_period, "\n")
  if (!is.null(x$followup_max)) cat("Follow-up times: ", x$followup_max, "\n")
  cat("Information times of looks: ", paste(format(x$timing, digits = 4), collapse = ", "), "\n")
  cat("Maximum information: ", paste(format(x$max_info, digits = 4), collapse = ", "), "\n")
  cat("Sample size group 1: ", x$n1, "\n")
  cat("Sample size group 2: ", x$n2, "\n")
  cat("Significance level: ", format(x$sig_level, digits = 4), "\n")
  cat("Power group sequential design: ", format(x$power_gs, digits = 4), "\n")
  cat("Power fixed design: ", format(x$power_fix, digits = 4), "\n")
  cat("\nProbabilities for rejection null hypothsis: \n")
  x$reject_prob <- as.data.frame(x$reject_prob)
  print(x$reject_prob, row.names = FALSE)
}

