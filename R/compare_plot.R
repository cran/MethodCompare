#' Plot used to visualize the recalibration of the new method after estimating
#' the bias
#' 
#' This function allows the visualization of the bias-corrected values (i.e.
#' recalibrated values, variable y1_corr) of the new measurement method.
#'
#' @param object list returned by \link{measure_compare} function.
#' 
#' @importFrom graphics title par points axis mtext box legend
#' 
#' @export
#'
#' @examples \donttest{
#' ### Load the data
#' data(data1)
#' ### Analysis
#' measure_model <- measure_compare(data1, nb_simul=100)
#' ### Plot the bias
#' compare_plot(measure_model)}
compare_plot <- function(object) {
  print("Generating Comparison Plot ...")
  data_sub <- object$data
  models <- object$models
  
  max <- max(data_sub$y2, data_sub$y1, data_sub$y1_corr, na.rm = TRUE)
  min <- min(data_sub$y2, data_sub$y1, data_sub$y1_corr, na.rm = TRUE)
  range <- max - min
  
  # Order data for plot
  data_sub <- data_sub[order(data_sub$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 2, 2) + 0.1)
  plot(data_sub$y2_hat, data_sub$y2, cex = 0.5, pch = 19,
       col = "darkgrey", axes = FALSE, xlab = "", ylab = "", 
       ylim = c(min - range * 0.1, max + range * 0.2))
  title(main = "Comparison of the methods", cex.main = 0.9)
  ### Add the y axis
  axis(2, col = "black", las = 1)
  mtext(sprintf("%s, %s_corr and %s", object$methods[1], object$methods[1], 
                object$methods[2]), side = 2, line = 2.5, cex = 0.8)
  box(col = "black")
  ### Add the x axis
  axis(1)
  mtext("True latent trait", side = 1, col = "black", line = 2, cex = 0.8)
  points(data_sub$y2_hat, data_sub$y1, pch = 19, col = "blue", cex = 0.5)
  points(data_sub$y2_hat, data_sub$y1_corr, pch = 19, col = "red",
         cex = 0.5)
  graphics::abline(c(0, 1), lwd = 2)
  graphics::abline(models[[4]]$coefficients, lwd = 2, lty = 1, col = "blue")
  graphics::abline(models[[6]]$coefficients, lwd = 2, lty = 2, col = "red")
  legend("topleft",
         legend = c(sprintf("Reference method (%s)", object$methods[2]),
                    sprintf("New method (%s)", object$methods[1]),
                    sprintf("Recalibrated new method (%s_corr)", 
                            object$methods[1])),
         pch = c(19, 19, 19), lty = c(1, 1, 2), col = c("black", "blue", "red"),
         y.intersp = 0.7, yjust = 0.2, bty = "n", cex = 0.8)
}