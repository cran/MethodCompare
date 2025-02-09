#' Estimation of the amount of bias of the new measurement method relative to
#' the reference method
#' 
#' `measure_compare()` implements the methodology reported in the paper:
#'  Taffé P. Effective plots to assess bias and precision in method comparison
#'  studies. Stat Methods Med Res 2018;27:1650-1660. Other relevant references:
#'  Taffé P, Peng M, Stagg V, Williamson T. Biasplot: A package to effective
#'  plots to assess bias and precision in method comparison studies. 
#'  Stata J 2017;17:208-221. Taffé P, Peng M, Stagg V, Williamson T.
#'  MethodCompare: An R package to assess bias and precision in method 
#'  comparison studies. Stat Methods Med Res 2019;28:2557-2565. 
#'  Taffé P, Halfon P, Halfon M. A new statistical methodology to assess bias
#'  and precision overcomes the defects of the Bland & Altman method.
#'  J Clin Epidemiol 2020;124:1-7. Taffé P. Assessing bias, precision, and
#'  agreement in method comparison studies. Stat Methods Med Res 2020;29:778-796. 
#'  Taffé P. When can the Bland-Altman limits of agreement method be used and
#'  when it should not be used. J Clin Epidemiol 2021; 137:176-181.
#'
#' @param data a required data frame containing the identification number of the 
#' subject (`id`), the measurement values from the new method (`y1`) and
#' those from the reference method (`y2`).
#' @param new an optional string. The column name containing the measurements of the new 
#' measurement method.
#' @param ref an optional string. The column name containing the measurements of the 
#' reference method (at least two measurements per subject).
#' @param id an optional string. The column name containing the subject 
#' identification numbers.
#' @param nb_simul an optional number. The number of simulations used for simultaneous
#' confidence bands.
#' @param if_value an optional number. Restrict the study to observed 
#' measurement greater than the provided value, i.e., `y1 >= if_value && y2 >= if_value`.
#'
#' @return The function returns a list with the following items:
#' * `models`: a list of models fitted in estimation procedure
#' * `data`: the original data frame with renamed columns and 
#'  additional computed data
#' * `sim_params`: estimated model coefficients used afterward
#' * `nb_simul`: the number of simulations used for confidence bands 
#'  simulations
#' * `bias`: differential and proportional biases for new method and the
#' associated 95 percent confidence intervals
#' * `methods`: a list of methods names provided by the user
#' 
#' @importFrom stats fitted
#' @importFrom stats vcov
#' @importFrom stats aggregate
#' @importFrom stats confint
#' @importFrom stats lm
#' 
#' @export
#'
#' @examples \donttest{
#' ### Load the data
#' data(data1)
#' ### Analysis
#' measure_model <- measure_compare(data1, nb_simul=100)}

measure_compare <- function(data, new = "y1", ref = "y2", id = "id", 
                            nb_simul = 1000, if_value = NULL) {
  print("Computing differential and proportional biases")
  print(paste("id variable:", id))
  print(paste("New method y variable:", new))
  print(paste("Reference method y variable:", ref))
  print(paste("Number of simulations set to", nb_simul))
  
  # Keep name as variable to return it
  y1 <- new
  y2 <- ref
  
  # Set variables as NULL for check() command
  y2_hat <- y1_corr <- sd_blup <- NULL
  
  # Data subset with only variables of interest
  data_sub <- data[, c(id, new, ref)]
  colnames(data_sub) <- c("id", "y1", "y2")
  
  if (!is.null(if_value)) {
    id_smaller <- data_sub[!(data_sub$y1 < if_value | data_sub$y2 < if_value), ]$id
    id_smaller <- unique(id_smaller)
    
    data_sub <- data_sub[data_sub$id %in% id_smaller, ]
  }
  
  data_sub$id <- factor(as.integer(data_sub$id))
  
  # Test if methods are equal and stop if they are 
  are_new_ref_equal <- identical(data_sub$y1, data_sub$y2)
  if(are_new_ref_equal) stop('Methods measurements are the same')
  
  # Drop rows where both y1 and y2 are missing
  data_sub <- data_sub[!(is.na(data_sub$y1) & is.na(data_sub$y2)), ]
  
  # Create a dataset for reference method by excluding missing value
  data_y2 <- data_sub[!is.na(data_sub$y2), ]
  
  # Model 1: Mixed model for BLUP estimate of x
  model_1 <- lme4::lmer(y2 ~ 1 + (1 | id), data = data_sub,
                        na.action = stats::na.exclude)
  data_sub$y2_hat <- fitted(model_1)
  
  # Computation of variances of blup
  
  v_blup <- var_blup(data_sub[!is.na(data_sub$y2), ], model = model_1)
  v_blup$v_blup <- v_blup$v_blup # + vcov(model_1)[1, 1]
  v_blup$sd_blup <- sqrt(v_blup$v_blup)
  
  data_sub <- merge(data_sub, v_blup, by = "id", all = TRUE)
  
  # Model 2: regression of y2 based on BLUP of x
  model_2 <- estimatr::lm_robust(y2 ~ 0 + y2_hat, data = data_sub,
                                 se_type = "stata", clusters = id)
  
  # Add residuals and absolute residuals of y2
  data_sub$resid_y2 <- data_sub$y2 - data_sub$y2_hat
  data_sub$resid_y2_abs <- abs(data_sub$resid_y2)
  
  # Model 3: Estimation of variance function for y2
  model_3 <- estimatr::lm_robust(resid_y2_abs ~ y2_hat, data = data_sub,
                                 se_type = "stata", clusters = id)
  
  # Model coefficients & variance-covariance matrix
  model_3_coef <- stats::coef(model_3)
  theta2_0e <- model_3_coef[1]
  theta2_1e <- model_3_coef[2]
  model_3_cov <- vcov(model_3)
  vtheta2_0e <- model_3_cov[1, 1]
  vtheta2_1e <- model_3_cov[2, 2]
  covtheta2e <- model_3_cov[1, 2]
  
  # Smooth standard deviation estimate and its variance
  fit_abs_res_y2 <- fitted(model_3)
  data_sub$fit_abs_res_y2 <- NA
  data_sub$fit_abs_res_y2[as.numeric(names(fit_abs_res_y2))] <- fit_abs_res_y2
  data_sub$sig_res_y2 <- data_sub$fit_abs_res_y2 * sqrt(pi / 2)
  data_sub$sig2_res_y2 <- data_sub$sig_res_y2^2
  data_sub$log_sig_res_y2 <- log(data_sub$sig_res_y2)
  
  v_fit_abs_res_y2 <- theta2_1e^2 * data_sub$v_blup + vtheta2_0e +
    vtheta2_1e * (data_sub$v_blup + data_sub$y2_hat^2) +
    2 * covtheta2e * data_sub$y2_hat
  
  v_sig_res_y2 <- (pi / 2) * v_fit_abs_res_y2
  data_sub$v_sig2_res_y2 <- pi^2 * data_sub$fit_abs_res_y2^2 * v_fit_abs_res_y2
  data_sub$se_sig_res_y2 <- sqrt(v_sig_res_y2)
  v_log_sig_res_y2 <- (1 / data_sub$fit_abs_res_y2^2) * v_fit_abs_res_y2
  data_sub$se_log_sig_res_y2 <- sqrt(v_log_sig_res_y2)
  
  # Model 4: Regression of y1 based on BLUP of x
  model_4 <- estimatr::lm_robust(y1 ~ y2_hat, data = data_sub, se_type = "stata",
                                 clusters = id)
  model_4_coef <- stats::coef(model_4)
  
  # Differential and proportional bias of new method
  bias <- cbind(model_4$coefficients, confint(model_4))
  data_sub$bias_y1 <- bias[1] + data_sub$y2_hat * (bias[2] - 1)
  rownames(bias) <- c("Differential bias", "Proportional bias")
  colnames(bias)[1] <- "Estimate"
  model_4_cov <- vcov(model_4)
  v_diff_bias <- model_4_cov[1, 1]
  v_prop_bias <- model_4_cov[2, 2]
  cov_bias <- model_4_cov[1, 2]
  
  # Bias variance
  data_sub$v_bias_y1 <- (bias[2] - 1)^2 * data_sub$v_blup + v_diff_bias +
    v_prop_bias * (data_sub$v_blup + data_sub$y2_hat^2) +
    2 * cov_bias * data_sub$y2_hat
  data_sub$se_bias_y1 <- sqrt(data_sub$v_bias_y1)
  
  # Compute residuals of model_4
  fit_model_4 <- fitted(model_4)
  data_sub$res_y1 <- NA
  data_sub$res_y1[as.numeric(names(fit_model_4))] <- fit_model_4
  data_sub$res_y1 <- data_sub$y1 - data_sub$res_y1
  data_sub$res_y1_abs <- abs(data_sub$res_y1)
  
  # Model 5: Variance estimation for y1 based on BLUP of x
  model_5 <- estimatr::lm_robust(res_y1_abs ~ y2_hat, data = data_sub,
                                 se_type = "stata", clusters = id)
  
  # Model coefficients & variance-covariance matrix
  model_5_coef <- stats::coef(model_5)
  theta1_0e <- model_5_coef[1]
  theta1_1e <- model_5_coef[2]
  model_5_cov <- vcov(model_5)
  vtheta1_0e <- model_5_cov[1, 1]
  vtheta1_1e <- model_5_cov[2, 2]
  covtheta1e <- model_5_cov[1, 2]
  
  # Smooth standard deviation estimate
  fit_abs_res_y1 <- fitted(model_5)
  data_sub$fit_abs_res_y1 <- NA
  data_sub$fit_abs_res_y1[as.numeric(names(fit_abs_res_y1))] <- fit_abs_res_y1

  data_sub$sig_res_y1 <- data_sub$fit_abs_res_y1 * sqrt(pi / 2)
  data_sub$sig2_res_y1 <- data_sub$sig_res_y1^2
  log_sig_res_y1 <- log(data_sub$sig_res_y1)
  
  data_sub$v_fit_abs_res_y1 <- theta1_1e^2 * data_sub$v_blup +
    vtheta1_0e +
    vtheta1_1e * (data_sub$v_blup + data_sub$y2_hat^2) +
    2 * covtheta1e * data_sub$y2_hat
  
  v_sig_res_y1 <- (pi / 2) * data_sub$v_fit_abs_res_y1
  data_sub$v_sig2_res_y1 <- pi^2 * data_sub$fit_abs_res_y1^2 *
    data_sub$v_fit_abs_res_y1
  v_log_sig_res_y1 <- (1 / (data_sub$fit_abs_res_y1^2)) *
    data_sub$v_fit_abs_res_y1
  
  # Corrected y1
  data_sub$y1_corr <- (data_sub$y1 - bias[1]) / bias[2]
  data_sub$bias <- bias[1] + data_sub$y2_hat * (bias[2] - 1)
  
  # Model 6: Regression on corrected y1 using BLUP of x
  model_6 <- lm(y1_corr ~ y2_hat, data = data_sub)

  fit_y1_corr <- fitted(model_6)
  data_sub$fit_y1_corr <- NA
  data_sub$fit_y1_corr[as.numeric(names(fit_y1_corr))] <- fit_y1_corr
  data_sub$res_y1_corr <- data_sub$y1_corr - data_sub$fit_y1_corr
  data_sub$res_y1_corr_abs <- abs(data_sub$res_y1_corr)
  
  # Model 7: Variance function estimation for corrected y1
  model_7 <- estimatr::lm_robust(res_y1_corr_abs ~ y2_hat, data = data_sub,
                                 se_type = "stata", clusters = id)
  model_7_coef <- stats::coef(model_7)
  
  # Model coefficients & variance-covariance matrix
  theta1_0_corre <- model_7_coef[1]
  theta1_1_corre <- model_7_coef[2]
  model_7_cov <- vcov(model_7)
  vtheta1_0_corre <- model_7_cov[1, 1]
  vtheta1_1_corre <- model_7_cov[2, 2]
  covtheta1_corre <- model_7_cov[1, 2]
  
  # Smooth standard estimate
  fit_res_y1_corr_abs <- fitted(model_7)
  data_sub$fit_res_y1_corr_abs <- NA
  data_sub$fit_res_y1_corr_abs[as.numeric(names(fit_res_y1_corr_abs))] <- fit_res_y1_corr_abs
  data_sub$sig_res_y1_corr <- data_sub$fit_res_y1_corr_abs * sqrt(pi / 2)
  data_sub$sig2_res_y1_corr <- data_sub$sig_res_y1_corr^2
  data_sub$log_sig_res_y1_corr <- log(data_sub$sig_res_y1_corr)
  
  v_fit_res_y1_corr_abs <- (theta1_1_corre^2) * data_sub$v_blup +
    vtheta1_0_corre + vtheta1_1_corre * (data_sub$v_blup +
                                           data_sub$y2_hat^2) +
    2 * covtheta1_corre * data_sub$y2_hat
  
  v_sig_res_y1_corr <- (pi / 2) * v_fit_res_y1_corr_abs
  data_sub$se_sig_res_y1_corr <- sqrt(v_sig_res_y1_corr)
  
  data_sub$v_sig2_res_y1_corr <- pi^2 *
    (data_sub$fit_res_y1_corr_abs^2) * v_fit_res_y1_corr_abs
  
  v_log_sig_res_y1_corr <- (1 / (data_sub$fit_res_y1_corr_abs^2)) *
    v_fit_res_y1_corr_abs
  
  data_sub$se_log_sig_res_y1_corr <- sqrt(v_log_sig_res_y1_corr)
  
  # Parameters for simulation
  sim_params <- list(
    model_3_coef = model_3_coef, model_3_cov = model_3_cov,
    model_4_coef = model_4_coef, model_4_cov = model_4_cov,
    model_5_coef = model_5_coef, model_5_cov = model_5_cov,
    model_7_coef = model_7_coef, model_7_cov = model_7_cov
  )
  
  return(list(models = list(model_1, model_2, model_3, model_4, model_5,
                            model_6,  model_7),
              data = data_sub[order(data_sub$id), ],
              sim_params = sim_params, nb_simul = nb_simul,
              bias = bias, methods = list(y1, y2)
  ))
}
