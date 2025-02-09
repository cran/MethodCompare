#' Plot the agreement after recalibration
#' 
#' This function draws the "agreement plot" after recalibration, which is used 
#' to visually appraise the degree of agreement between the new and reference 
#' methods, before recalibration of the new method.
#' It is obtained by graphing a scatter plot of `y1-y2` (difference of the methods) 
#' versus the BLUP of the latent trait, `x`, along with the bias and 95% limits
#' of agreement with their 95% simultaneous confidence bands. 
#' The function adds a second scale on the right axis, showing the percentage 
#' of agreement index.
#'
#' @inheritParams compare_plot
#' @param rarea if `TRUE`, draw the plot with shading areas between
#' the confidence bands.
#' 
#' @importFrom stats qnorm rnorm quantile
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
#' ### Plot the agreement after recalibration
#' agreement0(measure_model)}
agreement1 <- function(object, rarea = FALSE) {
  print("Generating Agreement Plot after recalibration ...")
  
  # Extract the objects from the output
  data_agg <- aggregate_data(object$data)
  data_new <- data_agg[!is.na(data_agg$y1), ]
  data_y1_y2 <- complete_data(object$data)
  params <- object$sim_params
  nb_simul <- object$nb_simul
  
  sig_d_corr <- sqrt(data_agg$sig2_res_y1_corr + data_agg$sig2_res_y2)
  sig2_d_corr <- sig_d_corr^2
  
  data_agg$pct_agreement_corr <- 1 - (qnorm(0.975) * sig_d_corr) /
    abs(data_agg$y2_hat)
  
  data_agg$loa_up_corr <- qnorm(0.975) * sig_d_corr
  data_agg$loa_lo_corr <- qnorm(0.025) * sig_d_corr
  
  cov_sig2_res_y2_sig2_res_y1_c <- pi^2 * data_agg$fit_abs_res_y2 *
    data_new$fit_res_y1_corr_abs * params$model_7_coef[2] *
    params$model_3_coef[2] * data_agg$v_blup
  
  v_sig_d_corr <- (1 / (4 * sig2_d_corr)) *
    (data_agg$v_sig2_res_y2 + data_new$v_sig2_res_y1_corr +
       2 * cov_sig2_res_y2_sig2_res_y1_c)
  
  v_loa_up_corr <- (qnorm(0.975)^2) * v_sig_d_corr
  se_loa_up_corr <- sqrt(v_loa_up_corr)
  
  v_loa_lo_corr <- (qnorm(0.975)^2) * v_sig_d_corr
  se_loa_lo_corr <- sqrt(v_loa_lo_corr)
  
  # Simulation parameters
  m1 <- matrix(params$model_7_coef, nrow = 1)
  v1 <- matrix(params$model_7_cov, ncol = 2)
  
  m2 <- matrix(params$model_3_coef, nrow = 1)
  v2 <- matrix(params$model_3_cov, ncol = 2)
  
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$y2_hat,
                      sd = data_agg$sd_blup)
    
    thetas1_corr_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m1, Sigma = v1)
    
    fit_abs_res_y1_corr_j <- thetas1_corr_j[, 1] + thetas1_corr_j[, 2] *
      blup_x_j
    
    sig_res_y1_corr_j <- fit_abs_res_y1_corr_j * sqrt(pi / 2)
    sig2_res_y1_corr_j <- sig_res_y1_corr_j^2
    
    v_fit_abs_res_y1_corr_j <- thetas1_corr_j[, 2]^2 * data_agg$v_blup +
      v1[1, 1] +
      v1[2, 2] * (data_agg$v_blup + blup_x_j^2) + 2 * v1[1, 2] * blup_x_j
    
    v_sig2_res_y1_corr_j <- pi^2 * fit_abs_res_y1_corr_j^2 *
      v_fit_abs_res_y1_corr_j
    
    thetas2_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m2, Sigma = v2)
    
    fit_abs_res_y2_j <- thetas2_j[, 1] + thetas2_j[, 2] * blup_x_j
    
    sig_res_y2_j <- fit_abs_res_y2_j * sqrt(pi / 2)
    sig2_res_y2_j <- sig_res_y2_j^2
    
    v_fit_abs_res_y2_j <- thetas2_j[, 2]^2 * data_agg$v_blup + v2[1, 1] +
      v2[2, 2] * (data_agg$v_blup + blup_x_j^2) + 2 * v2[1, 2] * blup_x_j
    
    v_sig2_res_y2_j <- pi^2 * fit_abs_res_y2_j^2 * v_fit_abs_res_y2_j
    
    sig2_d_corr_j <- sig2_res_y1_corr_j + sig2_res_y2_j
    sig_d_corr_j <- sqrt(sig2_d_corr_j)
    
    loa_up_corr_j <- qnorm(0.975) * sig_d_corr_j
    
    cov_s2_res_y2_s2_res_y1_c_j <- pi^2 * fit_abs_res_y2_j *
      fit_abs_res_y1_corr_j * thetas1_corr_j[, 2] * thetas2_j[, 2] *
      data_agg$v_blup
    
    v_sig_d_corr_j <- (1 / (4 * sig2_d_corr_j)) *
      (v_sig2_res_y2_j + v_sig2_res_y1_corr_j + 2 * cov_s2_res_y2_s2_res_y1_c_j)
    
    v_loa_up_corr_j <- qnorm(0.975)^2 * v_sig_d_corr_j
    
    d_j <- abs(loa_up_corr_j - data_agg$loa_up_corr) / sqrt(v_loa_up_corr_j)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value7 <- quantile(unlist(sim_max_d), c(0.95))
  
  data_agg$loa_up_corr_lo <- data_agg$loa_up_corr - crit_value7 *
    se_loa_up_corr
  data_agg$loa_up_corr_up <- data_agg$loa_up_corr + crit_value7 *
    se_loa_up_corr
    
  fp <- function(...) mfp::fp(...)
  
  frac_poly_loa_up_corr_lo <- mfp::mfp(loa_up_corr_lo ~ fp(y2_hat, df = 4),
                                  data = data_agg)
  frac_poly_loa_up_corr_up <- mfp::mfp(loa_up_corr_up ~ fp(y2_hat, df = 4),
                                  data = data_agg)
  
  data_agg$loa_up_corr_lo_fit <- predict(frac_poly_loa_up_corr_lo)
  data_agg$loa_up_corr_up_fit <- predict(frac_poly_loa_up_corr_up)
  
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$y2_hat,
                      sd = data_agg$sd_blup)
    
    thetas1_corr_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m1, Sigma = v1)
    
    fit_abs_res_y1_corr_j <- thetas1_corr_j[, 1] + thetas1_corr_j[, 2] *
      blup_x_j
    sig_res_y1_corr_j <- fit_abs_res_y1_corr_j * sqrt(pi / 2)
    sig2_res_y1_corr_j <- sig_res_y1_corr_j^2
    
    v_fit_abs_res_y1_corr_j <- thetas1_corr_j[, 2]^2 * data_agg$v_blup +
      v1[1, 1] +
      v1[2, 2] * (data_agg$v_blup + blup_x_j^2) + 2 * v1[1, 2] * blup_x_j
    
    v_sig2_res_y1_corr_j <- pi^2 * fit_abs_res_y1_corr_j^2 *
      v_fit_abs_res_y1_corr_j
    
    thetas2_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m2, Sigma = v2)
    
    fit_abs_res_y2_j <- thetas2_j[, 1] + thetas2_j[, 2] * blup_x_j
    sig_res_y2_j <- fit_abs_res_y2_j * sqrt(pi / 2)
    sig2_res_y2_j <- sig_res_y2_j^2
    
    v_fit_abs_res_y2_j <- thetas2_j[, 2]^2 * data_agg$v_blup + v2[1, 1] +
      v2[2, 2] * (data_agg$v_blup + blup_x_j^2) + 2 * v2[1, 2] * blup_x_j
    
    v_sig2_res_y2_j <- pi^2 * fit_abs_res_y2_j^2 * v_fit_abs_res_y2_j
    
    sig2_d_corr_j <- sig2_res_y1_corr_j + sig2_res_y2_j
    sig_d_corr_j <- sqrt(sig2_d_corr_j)
    
    loa_lo_corr_j <- - qnorm(0.975) * sig_d_corr_j
    
    cov_s2_res_y2_s2_res_y1_c_j <- pi^2 * fit_abs_res_y2_j *
      fit_abs_res_y1_corr_j * thetas1_corr_j[, 2] * thetas2_j[, 2] *
      data_agg$v_blup
    
    v_sig_d_corr_j <- (1 / (4 * sig2_d_corr_j)) *
      (v_sig2_res_y2_j + v_sig2_res_y1_corr_j + 2 * cov_s2_res_y2_s2_res_y1_c_j)
    
    v_loa_lo_corr_j <- qnorm(0.975)^2 * v_sig_d_corr_j
    
    d_j <- abs(loa_lo_corr_j - data_agg$loa_lo_corr) / sqrt(v_loa_lo_corr_j)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value8 <- quantile(unlist(sim_max_d), c(0.95))
  
  data_agg$loa_lo_corr_lo <- data_agg$loa_lo_corr - crit_value8 *
    se_loa_lo_corr
  data_agg$loa_lo_corr_up <- data_agg$loa_lo_corr + crit_value8 *
    se_loa_lo_corr
  
  frac_poly_loa_lo_corr_lo <- mfp::mfp(loa_lo_corr_lo ~ fp(y2_hat, df = 4),
                                  data = data_agg)
  frac_poly_loa_lo_corr_up <- mfp::mfp(loa_lo_corr_up ~ fp(y2_hat, df = 4),
                                  data = data_agg)
  
  data_agg$loa_lo_corr_lo_fit <- predict(frac_poly_loa_lo_corr_lo)
  data_agg$loa_lo_corr_up_fit <- predict(frac_poly_loa_lo_corr_up)
  
  # Compute min and max values for y-axis
  min_y <- min(data_agg$loa_up_corr_up_fit, data_agg$loa_up_corr_lo_fit,
               data_agg$loa_lo_corr_up_fit, data_agg$loa_lo_corr_lo_fit,
               na.rm = TRUE)
  max_y <- max(data_agg$loa_up_corr_up_fit, data_agg$loa_up_corr_lo_fit,
               data_agg$loa_lo_corr_up_fit, data_agg$loa_lo_corr_lo_fit,
               na.rm = TRUE)
  
  min_y <- floor(min_y)
  range <- max_y - min_y
  max_y <- ceiling(max_y + range * 0.2)
  
  data_y1_y2$diff <- data_y1_y2$y1_corr - data_y1_y2$y2
  
  # Order data for plot
  data_agg <- data_agg[order(data_agg$y2_hat), ]
  data_y1_y2 <- data_y1_y2[order(data_y1_y2$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the agreement after recalibration
  plot(data_y1_y2$y2_hat, data_y1_y2$diff, xlab = "",
       ylab = "", axes = FALSE, col = "darkgrey", ylim = c(min_y, max_y),
       cex = 0.5, pch = 19)
  title(main = "Agreement plot", cex.main = 0.9)
  
  # Bias
  graphics::abline(h = 0, col = "red", lty = 1)
  
  # Add the subtitle
  subtitle <- "(after recalibration)"
  mtext(subtitle, side = 3, cex = 0.8, line = .2)
  
  # 95% LoA
  points(data_agg$y2_hat, data_agg$loa_lo_corr, col = "dimgrey",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$loa_up_corr, col = "dimgrey",
         type = "l", lty = 2)
  
  # Lower confidence bands
  points(data_agg$y2_hat, data_agg$loa_lo_corr_up_fit, col = "orange",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$loa_lo_corr_lo_fit, col = "orange",
         type = "l", lty = 2)
  
  if (rarea) {
    polygon(
      c(data_agg$y2_hat, rev(data_agg$y2_hat)),
      c(data_agg$loa_lo_corr_lo_fit,
        rev(data_agg$loa_lo_corr_up_fit)),
      col = rgb(1, 0.5, 0, alpha = 0.2),
      border = NA
    )
  }
  
  # Upper confidence bands
  points(data_agg$y2_hat, data_agg$loa_up_corr_up_fit, col = "orange",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$loa_up_corr_lo_fit, col = "orange",
         type = "l", lty = 2)
  
  if (rarea) {
    polygon(
      c(data_agg$y2_hat, rev(data_agg$y2_hat)),
      c(data_agg$loa_up_corr_lo_fit,
        rev(data_agg$loa_up_corr_up_fit)),
      col = rgb(1, 0.5, 0, alpha = 0.2),
      border = NA
    )
  }
  
  # Left y-axis
  axis(2, col = "black", las = 1)
  mtext(sprintf("Difference (%s_corr - %s)", object$methods[1], 
                object$methods[2]), side = 2, line = 2.5, cex = 0.8)
  box(col = "black")
  
  # Add second plot: percentage agreement
  par(new = TRUE)
  plot(data_agg$y2_hat, data_agg$pct_agreement_corr, xlab = "",
       ylab = "", axes = FALSE, col = "blue", type = "l", lwd = 1,
       ylim = c(0, 1))
  
  # Right y-axis
  mtext("% of agreement", side = 4, col = "blue", line = 2.5, cex = 0.8)
  axis(4, col = "blue", col.axis = "blue", las = 1)
  
  # x-axis
  axis(1)
  mtext("True latent trait", side = 1, col = "black", line = 2, cex = 0.8)
  
  # Legend
  legend("top", legend = c("Bias", "95% LoA", "95% confidence limits",
                           "% of agreement"),
         pch = c(1, 19), col = c("red", "dimgrey", "orange", "blue"),
         pt.cex = c(0, 0), y.intersp = 0.7, yjust = 0.2, lty = c(1, 2, 2, 1),
         bty = "n", cex = 0.8)
}