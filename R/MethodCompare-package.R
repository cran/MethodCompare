#' Evaluating Bias and Precision in Method Comparison Studies
#' 
#' @description
#' The package "MethodCompare" allows one to assess bias, precision and agreement 
#' of a new measurement method with respect to a reference method (also called 
#' "reference standard"). It requires repeated measurements by at least one of 
#' the two measurement methods.
#' 
#' In this implementation, it is assumed by default that the reference method 
#' has repeated measurements and the new method may have as few as only one 
#' measurement per individual (The methodology can be adapted if you have more 
#' repeated measurements by the new method than by the reference method, see 
#' ref. below).
#' 
#' It implements the methodology developped in:
#' 
#' Taffé P. Effective plots to assess bias and precision in method comparison 
#' studies. Stat Methods Med Res 2018;27:1650-1660.
#' 
#' Taffé P. Assessing bias, precision, and agreement in method comparison 
#' studies. Stat Methods Med Res 2020;29:778-796.
#' 
#' For other relevant references:
#' 
#' Blomet T, Taffé P, MethodCompare: An extended suite of R commands to assess 
#' bias, precision, and agreement in method comparison studies.
#' To be published...
#' 
#' Taffé P, Peng M, Stagg V, Williamson T. Biasplot: A package to effective 
#' plots to assess bias and precision in method comparison studies. 
#' Stata J 2017;17:208-221.
#' 
#' Taffé P, Peng M, Stagg V, Williamson T. MethodCompare: An R package to 
#' assess bias and precision in method comparison studies. 
#' Stat Methods Med Res 2019;28:2557-2565.
#' 
#' Taffé P, Halfon P, Halfon M. A new statistical methodology to assess bias 
#' and precision overcomes the defects of the Bland & Altman method. J Clin Epidemiol 2020;124:1-7.
#' 
#' Taffé P. When can the Bland-Altman limits of agreement method be used and 
#' when it should not be used. J Clin Epidemiol 2021; 137:176-181.
#' 
#' Taffé P, Peng M, Stagg V, Williamson T. Extended biasplot command to assess 
#' bias, precision, and agreement in method comparison studies. 
#' Stata J 2023;23:97-118.
#' 
#' @details
#' The functions implemented in this package are the following:
#'  * \link{agreement0}: Plot the agreement before recalibration
#'  * \link{agreement1}: Plot the agreement after recalibration
#'  * \link{bias_plot}: Plot the bias and measurements
#'  * \link{compare_plot}: Plot used to visualize the recalibration of the new 
#'  method after estimating the bias
#'  * \link{measure_compare}: Estimation of the amount of bias of the new 
#'  measurement method relative to the reference method
#'  * \link{mse}: Plot the mean squared errors
#'  * \link{pct_agreement0}: Plot the percentage agreement before recalibration
#'  * \link{pct_agreement1}: Plot the percentage agreement after recalibration
#'  * \link{precision_plot}: Plot the precision of the methods
#'  * \link{sqrt_mse}: Plot the square root of the mean squared errors
#'  * \link{total_bias_plot}: Plot total bias
#' 
#' @keywords internal 
"_PACKAGE"