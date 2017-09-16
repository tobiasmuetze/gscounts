#' @name add_maxinfo
#' @title Maximum information
#' @description Add maximum information and critical values to the object gsnb
#' @keywords internal
#' @import stats
add_maxinfo <- function(x) {
  
  # Calculation of max information and critical values is split into three parts: 
  # (1) no futility
  # (2) non-binding futility
  # (3) binding futility
  
  k <- length(x$timing)
  log_effect <- log(x$rate1 / x$rate2) - log(x$ratio_H0)
  
  # (1) Maximum information for the case of no futility
  if (is.null(x$futility$type)) {
    ## Calculate critical values for efficacy testing 
    x$efficacy$critical[1] <- qnorm(x$efficacy$spend[1])
    for(i in 2:k) {
      x$efficacy$critical[i] <- cpp_calc_critical(r = 20, 
                                                  lower = x$efficacy$critical, 
                                                  upper = rep(Inf, times = i - 1), 
                                                  error_spend = x$efficacy$spend[i], 
                                                  information = x$timing[1:i], theta = 0, 
                                                  side = "lower")
    }
    ## Calculate maximum information for given power
    func <- function(max_info) { 
      reject_prob <- cpp_pmultinorm(r = 20, 
                                    lower = x$efficacy$critical,
                                    upper = rep(Inf, times = k),
                                    information = max_info * x$timing,
                                    theta = log_effect)
      1 - reject_prob - x$power
    }
    x$max_info <- uniroot(f = func, interval = c(0.1, 100), extendInt = "upX")$root
    return(x)
  }
  
  
  # (2) Maximum information for the case of non-binding futility
  if (x$futility$type == "nonbinding") {
    ## Calculate critical values for efficacy testing, independent of futility 
    x$efficacy$critical[1] <- qnorm(x$efficacy$spend[1])
    for(i in 2:k) {
      x$efficacy$critical[i] <- cpp_calc_critical(r = 20, 
                                                  lower = x$efficacy$critical, 
                                                  upper = rep(Inf, times = i - 1), 
                                                  error_spend = x$efficacy$spend[i], 
                                                  information = x$timing[1:i], theta = 0, 
                                                  side = "lower")
    }
    ## Calculate maximum information for given power
    func <- function(max_info) { 
      fut_critical <- numeric(k)
      ### Calculate critical values for non-binding futility
      fut_critical[1] <- qnorm(1 - x$futility$spend[1], 
                               mean = log_effect * sqrt(x$timing[1] * max_info))
      for(i in 2:k) {
        fut_critical[i] <- cpp_calc_critical(lower = x$efficacy$critical[1:(i-1)], 
                                             upper = fut_critical[1:(i-1)], 
                                             error_spend = x$futility$spend[i], 
                                             information = x$timing[1:i] * max_info, theta = log_effect, 
                                             side = "upper", r = 20)
        
      }
      x$efficacy$critical[k] - fut_critical[k]
    }
    x$max_info <- uniroot(f = func, interval = c(0.1, 100), extendInt = "upX")$root
    ## Calculate critical values for non-binding futility from maximum information
    x$futility$critical[1] <- qnorm(1 - x$futility$spend[1], 
                                    mean = log_effect * sqrt(x$timing[1] * x$max_info))
    for(i in 2:k) {
      x$futility$critical[i] <- cpp_calc_critical(r = 20, 
                                                  lower = x$efficacy$critical[1:(i-1)], 
                                      upper = x$futility$critical[1:(i-1)], 
                                      error_spend = x$futility$spend[i], 
                                      information = x$timing[1:i] * x$max_info, theta = log_effect, 
                                      side = "upper")
      
    }
    return(x)
  }
  
  # (3) Maximum information for the case of binding futility
  if (x$futility$type == "binding") {
    ## Calculate maximum information for given power
    func <- function(max_info) { 
      fut_critical <- eff_critical <- numeric(k)
      eff_critical[1] <- qnorm(x$efficacy$spend[1])
      fut_critical[1] <- qnorm(1 - x$futility$spend[1], 
                               mean = log_effect * sqrt(x$timing[1] * max_info))
      
      for(i in 2:k) {
        eff_critical[i] <- cpp_calc_critical(r = 20, 
                                             lower = eff_critical[1:(i-1)], 
                                 upper = fut_critical[1:(i-1)], 
                                 error_spend = x$efficacy$spend[i], 
                                 information = x$timing[1:i] * max_info, theta = 0, 
                                 side = "lower")
        fut_critical[i] <- cpp_calc_critical(r = 20, 
                                             lower = eff_critical[1:(i-1)], 
                                 upper = fut_critical[1:(i-1)], 
                                 error_spend = x$futility$spend[i], 
                                 information = x$timing[1:i] * max_info, theta = log_effect, 
                                 side = "upper")
      }
      eff_critical[k] - fut_critical[k]
    }
    x$max_info <- uniroot(f = func, interval = c(0.1, 100), extendInt = "upX")$root
    ## Calculate critical values for binding futility from maximum information
    x$futility$critical[1] <- qnorm(1 - x$futility$spend[1], 
                                    mean = log_effect * sqrt(x$timing[1] * x$max_info))
    x$efficacy$critical[1] <- qnorm(x$efficacy$spend[1])
    
    for(i in 2:k) {
      x$efficacy$critical[i] <- cpp_calc_critical(r = 20, 
                                                  lower = x$efficacy$critical[1:(i-1)], 
                                      upper = x$futility$critical[1:(i-1)], 
                                      error_spend = x$efficacy$spend[i], 
                                      information = x$timing[1:i] * x$max_info, theta = 0, 
                                      side = "lower")
      x$futility$critical[i] <- cpp_calc_critical(r = 20, 
                                                  lower = x$efficacy$critical[1:(i-1)], 
                                      upper = x$futility$critical[1:(i-1)], 
                                      error_spend = x$futility$spend[i], 
                                      information = x$timing[1:i] * x$max_info, theta = log_effect, 
                                      side = "upper")      
    }
    return(x)
  }
}