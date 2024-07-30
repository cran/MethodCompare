#' @importFrom stats rnorm
#' @importFrom stats quantile
precision_simulation <- function(object) {
  # Extract the objects from the output
  data_agg <- object$agg
  data_new <- object$new
  nb_simul <- object$nb_simul
  
  # Simulation for method 2
  # Retrieve useful params for simulation
  sim_params <- object$sim_params
  
  m2 <- matrix(data = sim_params$model_3_coef, nrow = 1)
  v2 <- matrix(data = sim_params$model_3_cov, ncol = 2)
  
  # Initiate variables
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$fitted_y2,
                             sd = data_agg$sd_blup)
    thetas_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m2, Sigma = v2)
    
    sig_res_y2_j <- (thetas_j[, 1] + thetas_j[, 2] * blup_x_j) * sqrt(pi / 2)
    
    v_fit_abs_res_y2_j <- thetas_j[, 2]^2 * data_agg$v_blup + v2[1, 1] +
      v2[2, 2] * (data_agg$v_blup + blup_x_j^2) + 2 * v2[1, 2] * blup_x_j
    v_sig_res_y2_j <- (pi / 2) * v_fit_abs_res_y2_j
    
    d_j <- abs(sig_res_y2_j - data_agg$sig_res_y2) / sqrt(v_sig_res_y2_j)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value1 <- quantile(unlist(sim_max_d), c(0.95))
  
  data_agg$sig_e2_lo <- exp(data_agg$log_sig_res_y2 - crit_value1 *
                              data_agg$se_log_sig_res_y2)
  data_agg$sig_e2_up <- exp(data_agg$log_sig_res_y2 + crit_value1 *
                              data_agg$se_log_sig_res_y2)
  
  fp <- function(...) mfp::fp(...)
  
  frac_poly_sig_e2_lo <- mfp::mfp(sig_e2_lo ~ fp(fitted_y2, df = 4), data = data_agg)
  frac_poly_sig_e2_up <- mfp::mfp(sig_e2_up ~ fp(fitted_y2, df = 4),
                                  data = data_agg)
  
  data_agg$sig_e2_lo_fit <- predict(frac_poly_sig_e2_lo)
  data_agg$sig_e2_up_fit <- predict(frac_poly_sig_e2_up)
  
  # Simulation for method 1
  # Retrieve useful params for simulation
  m1 <- matrix(data = sim_params$model_7_coef, nrow = 1)
  v1 <- matrix(data = sim_params$model_7_cov, ncol = 2)
  
  # Initiate variables
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$y2_hat,
                             sd = data_agg$sd_blup)
    thetas_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m1, Sigma = v1)
    
    sig_res_y1_corr_j <- (thetas_j[, 1] + thetas_j[, 2] * blup_x_j) *
      sqrt(pi / 2)
    
    v_fit_abs_res_y1_corr_j <- (thetas_j[, 2]^2) * data_agg$v_blup +
      v1[1, 1] + v1[2, 2] * (data_agg$v_blup + blup_x_j^2) +
      2 * v1[1, 2] * blup_x_j
    v_sig_res_y1_corr_j <- (pi / 2) * v_fit_abs_res_y1_corr_j
    
    d_j <- abs(sig_res_y1_corr_j -
                 data_agg$sig_res_y1_corr) / sqrt(v_sig_res_y1_corr_j)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value3 <- quantile(unlist(sim_max_d), c(0.95))
  
  data_agg$sig_e1_corr_lo <- exp(data_agg$log_sig_res_y1_corr -
                                   crit_value3 *
                                   data_agg$se_log_sig_res_y1_corr)
  data_agg$sig_e1_corr_up <- exp(data_agg$log_sig_res_y1_corr +
                                   crit_value3 *
                                   data_agg$se_log_sig_res_y1_corr)
  
  frac_poly_sig_e1_corr_lo <- mfp::mfp(sig_e1_corr_lo ~ fp(y2_hat, df = 4),
                                       data = data_agg)
  frac_poly_sig_e1_corr_up <- mfp::mfp(sig_e1_corr_up ~ fp(y2_hat, df = 4),
                                       data = data_agg)
  
  data_agg$sig_e1_corr_lo_fit <- predict(frac_poly_sig_e1_corr_lo)
  data_agg$sig_e1_corr_up_fit <- predict(frac_poly_sig_e1_corr_up)
  
  # Compute min and max values for y-axis
  min_y <- min(data_agg$sig_e2_lo_fit, data_agg$sig_e2_up_fit,
               data_agg$sig_e1_corr_lo_fit, data_agg$sig_e1_corr_up_fit,
               na.rm = TRUE)
  max_y <- max(data_agg$sig_e2_lo_fit, data_agg$sig_e2_up_fit,
               data_agg$sig_e1_corr_lo_fit, data_agg$sig_e1_corr_up_fit,
               na.rm = TRUE)
  
  min_y <- floor(min_y)
  range <- max_y - min_y
  max_y <- ceiling(max_y + range * 0.2)
  
  # Order data for plot
  data_agg <- data_agg[order(data_agg$y2_hat), ]
  data_new <- data_new[order(data_new$y2_hat), ]
  
  return(list(
    data_agg = data_agg, data_new = data_new,
    min_y = min_y, max_y = max_y
  ))
}