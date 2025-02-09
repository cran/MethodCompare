#' Plot the mean squared errors
#' 
#' This function draws the "MSE plot", which is used to compare the precision of 
#' the two measurement methods without recalibrating the new method. 
#' It is obtained by graphing the mean squared errors of `y1` (new method) and `y2` (reference 
#' method) versus the BLUP of the latent trait, `x`, along with their 95% 
#' simultaneous confidence bands.
#'
#' @inheritParams compare_plot
#' @param rarea if `TRUE`, draw the plot with shading areas between
#' the confidence bands.
#' 
#' @importFrom stats rnorm quantile
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
#' ### Plot the mean squared errors
#' mse(measure_model)}
mse <- function(object, rarea = FALSE) {
  print(" Generating MSE Plot ...")
  
  # Extract the objects from the output
  data_agg <- aggregate_data(object$data)
  params <- object$sim_params
  nb_simul <- object$nb_simul
  
  data_agg$mse1 <- data_agg$sig2_res_y1 + data_agg$bias_y1^2
  
  v_mse1 <- (pi^2) * (data_agg$fit_abs_res_y1^2) * data_agg$v_fit_abs_res_y1 +
    4 * (data_agg$bias_y1^2) * data_agg$v_bias_y1 +
    2 * pi * 2 * data_agg$fit_abs_res_y1 * data_agg$bias_y1 *
    params$model_5_coef[2] * (params$model_4_coef[2] - 1) *
    data_agg$v_blup
  se_mse1 <- sqrt(v_mse1)
  
  # Simulation parameters
  m1 <- matrix(params$model_5_coef, nrow = 1)
  v1 <- matrix(params$model_5_cov, ncol = 2)
  
  m2 <- matrix(params$model_3_coef, nrow = 1)
  v2 <- matrix(params$model_3_cov, ncol = 2)
  
  m3 <- matrix(params$model_4_coef, nrow = 1)
  v3 <- matrix(params$model_4_cov, ncol = 2)
  
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$y2_hat,
                      sd = data_agg$sd_blup)
    
    thetas1_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m1, Sigma = v1)
    
    fit_abs_res_y1_j <- thetas1_j[, 1] + thetas1_j[, 2] * blup_x_j
    
    sig_res_y1_j <- fit_abs_res_y1_j * sqrt(pi / 2)
    sig2_res_y1_j <- sig_res_y1_j^2
    
    v_fit_abs_res_y1_j <- (thetas1_j[, 2]^2) * data_agg$v_blup +
      v1[1, 1] + v1[2, 2] * (data_agg$v_blup + blup_x_j^2) +
      2 * v1[1, 2] * blup_x_j
    
    biases_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m3, Sigma = v3)
    
    v_bias_j <- (biases_j[, 2] - 1)^2 * data_agg$v_blup + v3[1, 1] +
      v3[2, 2] * (data_agg$v_blup + blup_x_j^2) + 2 * v3[1, 2] * blup_x_j
    
    bias_j <- biases_j[, 1] + blup_x_j * (biases_j[, 2] - 1)
    
    mse1_j <- sig2_res_y1_j + bias_j^2
    
    v_mse1_j <- (pi^2) * (fit_abs_res_y1_j^2) * v_fit_abs_res_y1_j +
      4 * (bias_j^2) * v_bias_j +
      2 * pi * fit_abs_res_y1_j * bias_j * thetas1_j[, 2] *
      (biases_j[, 2] - 1) * data_agg$v_blup
    
    d_j <- abs(mse1_j - data_agg$mse1) / sqrt(v_mse1_j)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value10 <- quantile(unlist(sim_max_d), c(0.95), na.rm = TRUE)
  
  data_agg$mse1_lo <- data_agg$mse1 - crit_value10 * se_mse1
  data_agg$mse1_up <- data_agg$mse1 + crit_value10 * se_mse1
  
  fp <- function(...) mfp::fp(...)
  
  frac_poly_mse1_lo <- mfp::mfp(mse1_lo ~ fp(y2_hat, df = 4), data = data_agg)
  frac_poly_mse1_up <- mfp::mfp(mse1_up ~ fp(y2_hat, df = 4), data = data_agg)
  
  data_agg$mse1_lo_fit <- predict(frac_poly_mse1_lo)
  data_agg$mse1_up_fit <- predict(frac_poly_mse1_up)
  
  data_agg$mse2 <- data_agg$sig2_res_y2
  
  v_mse2 <- data_agg$v_sig2_res_y2
  se_mse2 <- sqrt(v_mse2)
  
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$y2_hat,
                      sd = data_agg$sd_blup)
    
    thetas2_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m2, Sigma = v2)
    
    fit_abs_res_y2_j <- thetas2_j[, 1] + thetas2_j[, 2] * blup_x_j
    
    sig_res_y2_j <- fit_abs_res_y2_j * sqrt(pi / 2)
    sig2_res_y2_j <- sig_res_y2_j^2
    
    v_fit_abs_res_y2_j <- (thetas2_j[, 2]^2) * data_agg$v_blup +
      v2[1, 1] + v2[2, 2] * (data_agg$v_blup + blup_x_j^2) +
      2 * v2[1, 2] * blup_x_j
    
    v_sig2_res_y2_j <- pi^2 * fit_abs_res_y2_j^2 * v_fit_abs_res_y2_j
    
    mse2_j <- sig2_res_y2_j
    
    v_mse2_j <- v_sig2_res_y2_j
    
    d_j <- abs(mse2_j - data_agg$mse2) / sqrt(v_mse2_j)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value11 <- quantile(unlist(sim_max_d), c(0.95), na.rm = TRUE)
  
  data_agg$mse2_lo <- data_agg$mse2 - crit_value11 * se_mse2
  data_agg$mse2_up <- data_agg$mse2 + crit_value11 * se_mse2
  
  frac_poly_mse2_lo <- mfp::mfp(mse2_lo ~ fp(y2_hat, df = 4), data = data_agg)
  frac_poly_mse2_up <- mfp::mfp(mse2_up ~ fp(y2_hat, df = 4), data = data_agg)
  
  data_agg$mse2_lo_fit <- predict(frac_poly_mse2_lo)
  data_agg$mse2_up_fit <- predict(frac_poly_mse2_up)
  
  # Compute min and max values for y-axis
  min_y <- min(data_agg$mse1_lo_fit, data_agg$mse1_up_fit, data_agg$mse2_lo_fit, data_agg$mse2_up_fit, na.rm = TRUE)
  max_y <- max(data_agg$mse1_lo_fit, data_agg$mse1_up_fit, data_agg$mse2_lo_fit, data_agg$mse2_up_fit, na.rm = TRUE)
  
  min_y <- floor(min_y)
  range <- max_y - min_y
  max_y <- ceiling(max_y + range * 0.2)
  
  # Order data for plot
  data_agg <- data_agg[order(data_agg$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the MSE
  plot(data_agg$y2_hat, data_agg$mse1, xlab = "", ylab = "",
       axes = FALSE, col = "red", type = "l", lwd = 2, ylim = c(min_y, max_y))
  title(main = "MSE plot", cex.main = 0.9)
  
  # Confidence bands of MSE1
  points(data_agg$y2_hat, data_agg$mse1_lo_fit, col = "red",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$mse1_up_fit, col = "red",
         type = "l", lty = 2)
  
  if (rarea) {
    polygon(
      c(data_agg$y2_hat, rev(data_agg$y2_hat)),
      c(data_agg$mse1_lo_fit,
        rev(data_agg$mse1_up_fit)),
      col = rgb(1, 0, 0, alpha = 0.2),
      border = NA
    )
  }
  
  # MSE2
  points(data_agg$y2_hat, data_agg$mse2, col = "black",
         type = "l", lwd = 2)
  
  # Confidence bands of MSE2
  points(data_agg$y2_hat, data_agg$mse2_lo_fit, col = "black",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$mse2_up_fit, col = "black",
         type = "l", lty = 2)
  
  if (rarea) {
    polygon(
      c(data_agg$y2_hat, rev(data_agg$y2_hat)),
      c(data_agg$mse2_lo_fit,
        rev(data_agg$mse2_up_fit)),
      col = rgb(0, 0, 0, alpha = 0.2),
      border = NA
    )
  }
  
  # y-axis
  axis(2, col = "black", las = 1)
  mtext("MSE", side = 2, line = 2.5, cex = 0.8)
  box(col = "black")
  
  # x-axis
  axis(1)
  mtext("True latent trait", side = 1, col = "black", line = 2, cex = 0.8)
  
  # Legend
  legend("top",
         legend = c(sprintf("%s (Reference method)", object$methods[2]),
                    sprintf("%s (New method)", object$methods[1])),
         pch = c(1, 19), col = c("black", "red"), pt.cex = c(0, 0),
         y.intersp = 0.7, yjust = 0.2, lty = c(1, 1), bty = "n", cex = 0.8)
}