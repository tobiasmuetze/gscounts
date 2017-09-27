#' print.gsnb
#' @description print method for instance of class gsnb
#' @param x an object of class gsnb
#' @param ... optional arguments to print or plot methods
#' @export
print.gsnb <- function(x, ...) {
  title_text <- "\n\tGroup sequential trial with negative binomial outcomes"
  if (!is.null(x$futility$type)) 
    title_text <- paste0(title_text, " and with ", x$futility$type, " futility")
  cat(paste0(title_text, "\n\n"))
  cat("Distribution parameters\n")
  cat(" Rate group 1: ", x$rate1, "\n")
  cat(" Rate group 2: ", x$rate2, "\n")
  cat(" Dispersion parameter: ", x$dispersion, "\n")
  cat("Hypothesis testing\n")
  cat(" Rate ratio under null hypothesis: ", format(x$ratio_H0, digits = 4), "\n")
  cat(" Rate ratio under alternative: ", format(x$rate1 / x$rate2, digits = 4), "\n")
  cat(" Significance level: ", format(x$sig_level, digits = 4), "\n")
  cat(" Power group sequential design: ", format(x$power, digits = 4), "\n")
  cat(" Power fixed design: ", format(x$power_fix, digits = 4), "\n")
  cat(" Maximum information: ", paste(format(x$max_info, digits = 4), collapse = ", "), "\n")
  cat(" Number of looks: ", length(x$timing), "\n")
  cat(" Information times of looks: ", paste(format(x$timing, digits = 4), collapse = ", "), "\n")
  if (!is.null(x$calendar)) 
    cat(" Calendar times of looks: ", paste(format(x$calendar, digits = 4), collapse = ", "), "\n")
  cat("Sample size and study duration\n")
  cat(" Sample size group 1: ", x$n1, "\n")
  cat(" Sample size group 2: ", x$n2, "\n")
  if (!is.null(x$accrual_period)) cat(" Accrual period: ", x$accrual_period, "\n")
  if (!is.null(x$study_period)) cat(" Study duration: ", x$study_period, "\n")
  if (!is.null(x$followup_max)) cat(" Follow-up times: ", x$followup_max, "\n")
  
  # Output of critical values and spending
  critical_text <- "\nCritical values and spending at each data look"
  if (!is.null(x$futility)) {
    futility_text <- paste0(". The futility spending is ", x$futility$type, ".")
    critical_text <- paste0(critical_text, futility_text)
    tab_head <- sprintf("%55s %-s \n %7s %18s %12s %12s %12s %12s\n", "+------- Efficacy -------+", 
                        "+------- Futility -------+","Look", "Information time", "Spending", 
                        "Boundary", "Spending", "Boundary")
  } else {
    tab_head <- sprintf("%55s \n %7s %18s %12s %12s\n", "+------- Efficacy -------+",
                        "Look", "Information time", "Spending", 
                        "Boundary")
  }
  cat(paste0(critical_text, "\n"))
  cat(tab_head)
  for(i in seq_along(x$timing)) {
    row_string <- sprintf(" %7d %18.2g %12.5g %12.4f", i, x$timing[i], x$efficacy$spend[i], x$efficacy$critical[i])
    if (!is.null(x$futility)) 
      row_string <- paste0(row_string, 
                           sprintf(" %12.5g %12.4f", x$futility$spend[i], x$futility$critical[i]))
    cat(row_string)
    cat("\n")
  }
  
  # Probability of stopping for efficacy (i.e. rejection H0)
  cat("\nProbabilities of stopping for efficacy, i.e. for rejecting H0\n")
  print(x$stop_prob$efficacy, row.names = FALSE)
  
  
  # Probability of stopping for futility (i.e. accepting H0)
  if (!is.null(x$futility)) {
    cat("\nProbabilities of stopping for futility, i.e. for accepting H0\n")
    print(x$stop_prob$futility, row.names = FALSE)
  }  
  
  # Expected information
  cat("\nExpected information level\n")
  print(x$expected_info, row.names = FALSE)
  
  invisible(x)
}



#' print.nb
#' @description print method for instance of class nb
#' @param x an object of class nb
#' @param ... optional arguments to print or plot methods
#' @export
print.nb <- function(x, ...) {
  title_text <- "\n\tFixed sample design with negative binomial outcomes"
  cat(paste0(title_text, "\n\n"))
  cat("Distribution parameters\n")
  cat(" Rate group 1: ", x$rate1, "\n")
  cat(" Rate group 2: ", x$rate2, "\n")
  cat(" Dispersion parameter: ", x$dispersion, "\n")
  cat("Hypothesis testing\n")
  cat(" Rate ratio under null hypothesis: ", format(x$ratio_H0, digits = 4), "\n")
  cat(" Rate ratio under alternative: ", format(x$rate1 / x$rate2, digits = 4), "\n")
  cat(" Significance level: ", format(x$sig_level, digits = 4), "\n")
  cat(" Power fixed design: ", format(x$power, digits = 4), "\n")
  cat(" Maximum information: ", format(x$max_info, digits = 4), "\n")
  cat("Sample size and study duration\n")
  cat(" Sample size group 1: ", x$n1, "\n")
  cat(" Sample size group 2: ", x$n2, "\n")
  if (!is.null(x$accrual_period)) cat(" Accrual period: ", x$accrual_period, "\n")
  if (!is.null(x$study_period)) cat(" Study duration: ", x$study_period, "\n")
  if (!is.null(x$followup_max)) cat(" Follow-up times: ", x$followup_max, "\n")
  
  invisible(x)
}

