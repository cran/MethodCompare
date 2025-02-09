#' Plot the precision of the methods
#' 
#' This function draws the "precision plot", which allows the visual comparison 
#' of the precision (i.e. standard deviation) of the new measurement method with 
#' the reference standard by creating a scatter plot of the estimated standard 
#' deviations, along with their 95% simultaneous confidence bands, against the 
#' best linear prediction (BLUP) of the true latent trait, `x`.
#'
#' @inheritParams compare_plot
#' @param object2 (optional) returned by \link{measure_compare} function. 
#' If provided, will plot a second precision estimate.
#' @param log if `TRUE`, guarantee the simultaneous confidence bands around the 
#' standard deviation of measurement errors to be strictly positive.
#' @param rarea if `TRUE`, draw the plot with shading areas between
#' the confidence bands.
#' 
#' @importFrom graphics title par points axis mtext box legend polygon
#' @importFrom grDevices rgb
#' 
#' @export
#'
#' @examples \donttest{
#' ### Load the data
#' data(data1)
#' ### Analysis
#' measure_model <- measure_compare(data1, nb_simul=100)
#' ### Plot the precision of the two methods
#' precision_plot(measure_model)}
precision_plot <- function(object, object2 = NULL, log = FALSE, rarea = FALSE) {
  print("Generating Precision Plot ...")
  
  two_objects <- !is.null(object2)
  
  # Simulate for the first dataset
  object_sim <- precision_simulation(object, log = log)
  min_y <- object_sim$min_y
  max_y <- object_sim$max_y
  
  # Simulate for the second dataset (if there is one)
  if (two_objects) {
    object2_sim <- precision_simulation(object2, log = log)
    min_y <- min(min_y, object2_sim$min_y)
    max_y <- max(max_y, object2_sim$max_y)
  }
  
  par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the precision
  plot(object_sim$data_agg$y2_hat, object_sim$data_agg$sig_res_y1_corr,
       xlab = "", ylab = "", axes = FALSE, col = "red", type = "l", lwd = 2,
       ylim = c(min_y, max_y))
  title(main = "Precision plot", cex.main = 0.9)
  
  if (two_objects) {
    points(object2_sim$data_agg$y2_hat,
           object2_sim$data_agg$sig_res_y1_corr, xlab = "", ylab = "",
           col = "blue", type = "l", lwd = 2)
  }
  
  # Add the subtitle
  subtitle <- "with 95% confidence bands (after recalibration)"
  mtext(subtitle, side = 3, cex = 0.8, line = .2)
  
  # Confidence bands
  points(object_sim$data_agg$y2_hat,
         object_sim$data_agg$sig_e1_corr_lo_fit, col = "red",
         type = "l", lty = 2)
  points(object_sim$data_agg$y2_hat,
         object_sim$data_agg$sig_e1_corr_up_fit, col = "red",
         type = "l", lty = 2)
  
  if (rarea) {
    polygon(
      c(object_sim$data_agg$y2_hat, rev(object_sim$data_agg$y2_hat)),
      c(object_sim$data_agg$sig_e1_corr_lo_fit,
        rev(object_sim$data_agg$sig_e1_corr_up_fit)),
      col = rgb(1, 0, 0, alpha = 0.2),
      border = NA
    )
  }
  
  if (two_objects) {
    points(object2_sim$data_agg$y2_hat,
           object2_sim$data_agg$sig_e1_corr_lo_fit, col = "blue",
           type = "l", lty = 2)
    points(object2_sim$data_agg$y2_hat,
           object2_sim$data_agg$sig_e1_corr_up_fit, col = "blue",
           type = "l", lty = 2)
    
    if (rarea) {
      polygon(
        c(object2_sim$data_agg$y2_hat, rev(object2_sim$data_agg$y2_hat)),
        c(object2_sim$data_agg$sig_e1_corr_lo_fit,
          rev(object2_sim$data_agg$sig_e1_corr_up_fit)),
        col = rgb(0, 0, 1, alpha = 0.2),
        border = NA
      )
    }
  }
  
  # y2
  points(object_sim$data_agg$y2_hat, 
         object_sim$data_agg$sig_res_y2, col = "black",
         type = "l", lwd = 2)
  
  # Confidence bands
  points(object_sim$data_agg$y2_hat,
         object_sim$data_agg$sig_e2_lo_fit, col = "black",
         type = "l", lty = 2)
  points(object_sim$data_agg$y2_hat,
         object_sim$data_agg$sig_e2_up_fit, col = "black",
         type = "l", lty = 2)
  
  if (rarea) {
    polygon(
      c(object_sim$data_agg$y2_hat, rev(object_sim$data_agg$y2_hat)),
      c(object_sim$data_agg$sig_e2_lo_fit,
        rev(object_sim$data_agg$sig_e2_up_fit)),
      col = rgb(0, 0, 0, alpha = 0.2),
      border = NA
    )
  }
  
  # y-axis
  axis(2, col = "black", las = 1)
  mtext("Standard deviation of the measurement errors", side = 2, line = 2.5, cex = 0.8)
  box(col = "black")
  
  # x-axis
  axis(1)
  mtext("True latent trait", side = 1, col = "black", line = 2, cex = 0.8)
  
  # Legend
  if (!two_objects) {
    legend("top", legend = c(sprintf("%s (Reference method)", object$methods[2]),
                             sprintf("%s_corr (Recalibrated new method)",
                                     object$methods[1])),
           pch = c(1, 19), col = c("black", "red"), pt.cex = c(0, 0),
           y.intersp = 0.7, yjust = 0.2, lty = c(1, 1), bty = "n",
           cex = 0.8)
  } else {
    legend("top", legend = c(sprintf("%s (Reference method)", object$methods[2]),
                             sprintf("%s_corr (Recalibrated new method: %s)",
                                     object$methods[1],
                                     deparse(substitute(object))),
                             sprintf("%s_corr (Recalibrated new method: %s)",
                                     object2$methods[1],
                                     deparse(substitute(object2)))),
           pch = c(1, 19), col = c("black", "red", "blue"), pt.cex = c(0, 0),
           y.intersp = 0.7, yjust = 0.2, lty = c(1, 1, 1), bty = "n", cex = 0.8)
  }
}