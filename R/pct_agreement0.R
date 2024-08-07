#' Plot the percentage agreement before recalibration
#' 
#' This function draws the "percentage agreement plot" before recalibration, 
#' which shows the amount of percentage agreement.
#' It is obtained by graphing the percentage agreement index before recalibration 
#' versus the BLUP of the latent trait, `x`, along with its 95% simultaneous 
#' confidence bands.
#'
#' @inheritParams compare_plot
#' 
#' @importFrom stats qnorm rnorm quantile var
#' @importFrom graphics title par points axis mtext box legend
#' 
#' @export
#'
#' @examples \donttest{
#' ### Load the data
#' data(data1)
#' ### Analysis
#' measure_model <- measure_compare(data1, nb_simul=100)
#' ### Plot the percentage agreement without recalibration
#' pct_agreement0(measure_model)}
pct_agreement0 <- function(object) {
  print("Generating Percentage Agreement Plot without recalibration ...")
  
  # Extract the objects from the output
  data_old <- object$ref
  data_agg <- object$agg
  params <- object$sim_params
  nb_simul <- object$nb_simul
  
  nb_ind <- unique(data_old$id)
  
  sig2_d <- data_agg$sig2_res_y1 + data_agg$sig2_res_y2
  sig_d <- sqrt(sig2_d)
  
  data_agg$pct_agreement <- 1 - (qnorm(0.975) * sig_d + abs(data_agg$bias_y1)) /
    abs(data_agg$fitted_y2)
  
  # Simulation parameters
  m1 <- matrix(params$model_5_coef, nrow = 1)
  v1 <- matrix(params$model_5_cov, ncol = 2)
  
  m2 <- matrix(params$model_3_coef, nrow = 1)
  v2 <- matrix(params$model_3_cov, ncol = 2)
  
  m3 <- matrix(params$model_4_coef, nrow = 1)
  v3 <- matrix(params$model_4_cov, ncol = 2)
  
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in nb_ind) {
    m_blup_x_j <- min(data_old[data_old$id == j, ]$fitted_y2)
    v_blup_x_j <- min(data_old[data_old$id == j, ]$v_blup)
    sd_blup_x_j <- sqrt(v_blup_x_j)
    
    blup_x_j <- rnorm(10000, mean = m_blup_x_j,
                      sd = sd_blup_x_j)
    
    thetas1_j <- rockchalk::mvrnorm(10000, mu = m1, Sigma = v1)
    thetas2_j <- rockchalk::mvrnorm(10000, mu = m2, Sigma = v2)
    biases_j <- rockchalk::mvrnorm(10000, mu = m3, Sigma = v3)
    
    sig_d_j <- sqrt((pi / 2) * (thetas1_j[, 1] + thetas1_j[, 2] * blup_x_j)^2 +
                      (pi / 2) * (thetas2_j[, 1] + thetas2_j[, 2] * blup_x_j)^2)
    bias_j <- biases_j[, 1] + blup_x_j * (biases_j[, 2] - 1)
    
    pct_agreement_j <- 1 - (qnorm(0.975) * sig_d_j + abs(bias_j)) /
      abs(blup_x_j)
    
    data_agg$v_pct_agreement[data_agg$id == j] <- var(pct_agreement_j)
  }
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$fitted_y2,
                      sd = data_agg$sd_blup)
    
    thetas1_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m1, Sigma = v1)
    
    sig_res_y1_j <- (thetas1_j[, 1] + thetas1_j[, 2] * blup_x_j) * sqrt(pi / 2)
    sig2_res_y1_j <- sig_res_y1_j^2
    
    thetas2_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m2, Sigma = v2)
    
    sig_res_y2_j <- (thetas2_j[, 1] + thetas2_j[, 2] * blup_x_j) * sqrt(pi / 2)
    sig2_res_y2_j <- sig_res_y2_j^2
    
    sig_d_j <- sqrt(sig2_res_y1_j + sig2_res_y2_j)
    
    biases_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m3, Sigma = v3)
    
    bias_j <- biases_j[, 1] + blup_x_j * (biases_j[, 2] - 1)
    
    pct_agreement_j <- 1 - (qnorm(0.975) * sig_d_j + abs(bias_j)) /
      abs(blup_x_j)
    
    d_j <- abs(pct_agreement_j - data_agg$pct_agreement) /
      sqrt(data_agg$v_pct_agreement)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value6 <- quantile(unlist(sim_max_d), c(0.95))
  
  data_agg$pct_agreement_lo <- data_agg$pct_agreement - crit_value6 *
    sqrt(data_agg$v_pct_agreement)
  data_agg$pct_agreement_up <- data_agg$pct_agreement + crit_value6 *
    sqrt(data_agg$v_pct_agreement)
  
  fp <- function(...) mfp::fp(...)
  
  frac_poly_pct_agreement_lo <- mfp::mfp(pct_agreement_lo ~ fp(fitted_y2, df = 4),
                                    data = data_agg)
  frac_poly_pct_agreement_up <- mfp::mfp(pct_agreement_up ~ fp(fitted_y2, df = 4),
                                    data = data_agg)
  
  data_agg$pct_agreement_lo_fit <- predict(frac_poly_pct_agreement_lo)
  data_agg$pct_agreement_up_fit <- predict(frac_poly_pct_agreement_up)
  
  # Compute min and max values for y-axis
  min_y <- min(data_agg$pct_agreement_lo_fit, data_agg$pct_agreement_up_fit,
               data_agg$pct_agreement, na.rm = TRUE)
  max_y <- max(data_agg$pct_agreement_lo_fit, data_agg$pct_agreement_up_fit,
               data_agg$pct_agreement, na.rm = TRUE)
  
  range <- max_y - min_y
  max_y <- max_y + range * 0.2
  
  # Order data for plot
  data_agg <- data_agg[order(data_agg$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the percentage agreement with no recalibration
  plot(data_agg$y2_hat, data_agg$pct_agreement, xlab = "", ylab = "",
       axes = FALSE, col = "blue", ylim = c(max(min_y, 0), min(max_y, 1)),
       type = "l")
  title(main = "Percentage agreement plot", cex.main = 0.9)
  
  # Add the subtitle
  subtitle <- "(no recalibration)"
  mtext(subtitle, side = 3, cex = 0.8)
  
  # 95% confidence bands
  points(data_agg$y2_hat, data_agg$pct_agreement_lo_fit, col = "blue",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$pct_agreement_up_fit, col = "blue",
         type = "l", lty = 2)
  
  # y-axis
  axis(2, col = "black", las = 1)
  mtext("Percentage agreement", side = 2, line = 2)
  box(col = "black")
  
  # x-axis
  axis(1)
  mtext("BLUP of x", side = 1, col = "black", line = 2)
  
  # Legend
  legend("top", legend = c("% agreement", "95% CB"),
         pch = c(1, 19), col = c("blue", "blue"),
         pt.cex = c(0, 0), y.intersp = 0.7, yjust = 0.2, lty = c(1, 2),
         bty = "n", cex = 0.8)
}