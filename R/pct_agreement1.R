#' Plot the percentage agreement after recalibration
#' 
#' This function draws the "percentage agreement plot" after recalibration, 
#' which shows the amount of percentage agreement.
#' It is obtained by graphing the percentage agreement index after recalibration 
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
#' ### Plot the percentage agreement after recalibration
#' pct_agreement0(measure_model)}
pct_agreement1 <- function(object) {
  print("Generating Percentage Agreement Plot after recalibration ...")
  
  # Extract the objects from the output
  data_old <- object$data
  data_agg <- aggregate_data(object$data)
  params <- object$sim_params
  nb_simul <- object$nb_simul
  
  nb_ind <- unique(data_old$id)
  
  sig_d_corr <- sqrt(data_agg$sig2_res_y1_corr + data_agg$sig2_res_y2)
  sig2_d_corr <- sig_d_corr^2
  
  data_agg$pct_agreement_corr <- 1 - (qnorm(0.975) * sig_d_corr) /
    abs(data_agg$y2_hat)
  
  # Simulation parameters
  m1 <- matrix(params$model_7_coef, nrow = 1)
  v1 <- matrix(params$model_7_cov, ncol = 2)
  
  m2 <- matrix(params$model_3_coef, nrow = 1)
  v2 <- matrix(params$model_3_cov, ncol = 2)
  
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in nb_ind) {
    m_blup_x_j <- min(data_old[data_old$id == j, ]$y2_hat, na.rm = TRUE)
    v_blup_x_j <- min(data_old[data_old$id == j, ]$v_blup, na.rm = TRUE)
    sd_blup_x_j <- sqrt(v_blup_x_j)
    
    blup_x_j <- abs(rnorm(1000, mean = m_blup_x_j, sd = sd_blup_x_j))
    
    thetas1_corr_j <- rockchalk::mvrnorm(1000, mu = m1, Sigma = v1)
    thetas2_j <- rockchalk::mvrnorm(1000, mu = m2, Sigma = v2)
    
    sig_d_corr_j <- sqrt((pi / 2) * (thetas1_corr_j[, 1] + thetas1_corr_j[, 2] *
                                       blup_x_j)^2 +
                           (pi / 2) * (thetas2_j[, 1] + thetas2_j[, 2] *
                                         blup_x_j)^2)
    
    pct_agreement_corr_j <- 1 - (qnorm(0.975) * sig_d_corr_j) / blup_x_j
    pct_agreement_corr_j <- pmax(pct_agreement_corr_j, 0)
    
    data_agg$v_pct_agreement_corr[data_agg$id == j] <- var(
      pct_agreement_corr_j
    )
  }
  
  for (j in 1:nb_simul) {
    blup_x_j <- abs(rnorm(dim(data_agg)[1], mean = data_agg$y2_hat,
                      sd = data_agg$sd_blup))
    
    thetas1_corr_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m1, Sigma = v1)
    thetas2_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m2, Sigma = v2)
    
    sig_d_corr_j <- sqrt((pi / 2) * (thetas1_corr_j[, 1] + thetas1_corr_j[, 2] *
                                       blup_x_j)^2 +
                           (pi / 2) * (thetas2_j[, 1] + thetas2_j[, 2] *
                                         blup_x_j)^2)
    
    pct_agreement_corr_j <- 1 - (qnorm(0.975) * sig_d_corr_j) / blup_x_j
    pct_agreement_corr_j <- pmax(pct_agreement_corr_j, 0)
    
    d_j <- abs(pct_agreement_corr_j - data_agg$pct_agreement_corr) /
      sqrt(data_agg$v_pct_agreement_corr)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value9 <- quantile(unlist(sim_max_d), c(0.95))
  
  data_agg$pct_agreement_corr_lo <- data_agg$pct_agreement_corr -
    crit_value9 * sqrt(data_agg$v_pct_agreement_corr)
  data_agg$pct_agreement_corr_up <- data_agg$pct_agreement_corr +
    crit_value9 * sqrt(data_agg$v_pct_agreement_corr)
  
  fp <- function(...) mfp::fp(...)
  
  frac_poly_pct_agreement_c_lo <- mfp::mfp(pct_agreement_corr_lo ~ fp(y2_hat,
                                                                 df = 4),
                                      data = data_agg)
  frac_poly_pct_agreement_c_up <- mfp::mfp(pct_agreement_corr_up ~ fp(y2_hat,
                                                                 df = 4),
                                      data = data_agg)
  
  data_agg$pct_agreement_c_lo_fit <- predict(frac_poly_pct_agreement_c_lo)
  data_agg$pct_agreement_c_up_fit <- predict(frac_poly_pct_agreement_c_up)
  
  # Compute min and max values for y-axis
  min_y <- min(data_agg$pct_agreement_c_lo_fit, data_agg$pct_agreement_c_up_fit,
               data_agg$pct_agreement_corr, na.rm = TRUE)
  max_y <- max(data_agg$pct_agreement_c_lo_fit, data_agg$pct_agreement_c_up_fit,
               data_agg$pct_agreement_corr, na.rm = TRUE)
  
  range <- max_y - min_y
  max_y <- max_y + range * 0.2
  
  # Order data for plot
  data_agg <- data_agg[order(data_agg$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the percentage agreement after recalibration
  plot(data_agg$y2_hat, data_agg$pct_agreement_corr, xlab = "",
       ylab = "", axes = FALSE, col = "blue",
       ylim = c(min_y, max_y), type = "l")
  title(main = "Percentage agreement plot", cex.main = 0.9)
  
  # Add the subtitle
  subtitle <- "(after recalibration)"
  mtext(subtitle, side = 3, cex = 0.8, line = .2)
  
  # 95% confidence bands
  points(data_agg$y2_hat, data_agg$pct_agreement_c_lo_fit, col = "blue",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$pct_agreement_c_up_fit, col = "blue",
         type = "l", lty = 2)
  
  # y-axis
  axis(2, col = "black", las = 1)
  mtext("Percentage agreement", side = 2, line = 2.75, cex = 0.8)
  box(col = "black")
  
  # x-axis
  axis(1)
  mtext("True latent trait", side = 1, col = "black", line = 2, cex = 0.8)
  
  # Legend
  legend("top", legend = c("% agreement", "95% CB"),
         pch = c(1, 19), col = c("blue", "blue"),
         pt.cex = c(0, 0), y.intersp = 0.7, yjust = 0.2, lty = c(1, 2),
         bty = "n", cex = 0.8)
}