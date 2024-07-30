#' Simulated dataset 3
#' 
#' In the simulated dataset 3, each subject has 10 to 20 measurement values
#' from the new method and 10 to 20 measurement values from the reference method.
#' Compared to the reference method, the new method has differential bias of 1 
#' and proportional bias of 0.9. Variance of the new method is smaller than that
#' for the reference method.
#' 
#' @format ## `data3`
#' An object of class `data.frame` with 1682 rows and 3 columns
#' 
#' @usage data3
#' 
#' @details
#' A data frame with 3 variables:
#' \describe{
#'  \item{`id`}{identification number for subjects}
#'  \item{`y1`}{values from the new measuremment method}
#'  \item{`y2`}{values from the reference method}
#' }
#' Dataset 1 was created based on the following equations:
#' \deqn{y_{1i}=1+0.9x_i+\varepsilon_{1i},\quad \varepsilon_{1i} \mid x_i \sim 
#' N(0,(1+0.04x_i)^2)}{y_1i=4+0.8*x_i+\epsilon_1i,	\epsilon_1i | x_i ~ 
#' N(0,(1+0.04*x_i)^2)}
#' \deqn{y_{2i}=x_i+\varepsilon_{2i},\quad \varepsilon_{2i} \mid x_i \sim 
#' N(0,(1.75+0.08x_i)^2)}{y_2i=x_i+\epsilon_2i,	\epsilon_2i | x_i ~ 
#' N(0,(1.75+0.08*x_i)^2)}
#' \deqn{x_i\sim Uniform[20-100]}
#'
#' for \eqn{i=1,\ldots,100} and the number of repeated measurements for each 
#' subject \eqn{i} from the reference standard was \eqn{n_{2i} \sim Uniform[10,20]}
#' and \eqn{n_{1i} \sim Uniform[10,20]} for the new measurement method.
"data3"