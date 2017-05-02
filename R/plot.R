#' plot.GSNBpower
#' @description Plots an object of class GSNBdesign
#' @details The following entries of the objects of class \emph{GSNBdesign}
#' are considered to calculate the power: \code{ratio_H0}, \code{ratio_H1},
#' \code{timing}, \code{sig_level}, \code{esf}
#' @param x an object of class \emph{GSNBpower}
#' @return Returns a gg-plot, i.e. an instance with classes \emph{gg} and \emph{ggplot}
#' @import ggplot2
plot.GSNBpower <- function(x) {

  ratio_H0 <- x$ratio_H0
  ratio_H1 <- x$ratio_H1
  timing <- x$timing
  sig_level <- x$sig_level
  esf <- x$esf

  # Calculate the information corresponding to the power limits of 2*siglevel und 0.99
  max_info_min <- power_gsnb(ratio_H1 = ratio_H1, power_gs = 2*sig_level, timing = timing,
                             esf = esf, ratio_H0 = ratio_H0, sig_level = sig_level)$max_info
  max_info_max <- power_gsnb(ratio_H1 = ratio_H1, power_gs = 0.99, timing = timing,
                             esf = esf, ratio_H0 = ratio_H0, sig_level = sig_level)$max_info

  # Function for power calculation
  body_powerGS <- quote({
    power_gsnb(ratio_H1 = ratio_H1, max_info = max_info, timing = timing,
               esf = esf, ratio_H0 = ratio_H0, sig_level = sig_level)$power_gs
  })

  # Define maxinfo for x-axis and calculate power
  maxinfo_plot <- seq(max_info_min, max_info_max, length.out = 100)
  powergs_plot <- sapply(X = maxinfo_plot, FUN = function(max_info) {eval(body_powerGS)})


  ggplot(mapping = aes(x = maxinfo_plot, y = powergs_plot)) +
    geom_line() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_x_continuous(breaks = round(seq(min(maxinfo_plot), max(maxinfo_plot), length.out = 10))) +
    ylab(label = "Power group sequential design") +
    xlab(expression(paste("Maximum information level ", I[max]))) +
    theme_bw() +
    theme(legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          text = element_text(size = 13))
}
