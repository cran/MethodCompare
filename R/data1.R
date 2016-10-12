#' Simulated dataset 1
#'
#' In the simulated dataset 1, each subject has 1 measurement value from the new method
#' and 10 to 15 measurement values from the reference method.
#' Compared to the reference method, the new method has
#' differential bias of -4 and proportional bias of 1.2. Variance of the new method
#' is smaller than that for the referene method.
#'
#'  @format A data frame with three variables:
#' \describe{
#' \item{\code{id}}{identification number for subjects}
#' \item{\code{y1}}{values from the new measurement method}
#' \item{\code{y2}}{values from the reference measurement method}
#' }
#' @details
#' Dataset 1 was created based on the following equations:
#' \deqn{y_{1i}=-4+1.2x_i+\varepsilon_{1i}, \varepsilon_{1i} \mid x_i \sim N(0,(1+0.1x_i)^2)}{
#' y_1i=-4+1.2*x_i+\epsilon_1i,\epsilon_1i|x_i~N(0,1+0.1*x_i)^2)}
#' \deqn{y_{2ij}=x_i+\varepsilon_{2ij},\varepsilon_{2ij} \mid x_i\sim N(0,(2+0.2x_i)^2)}{
#' y_2ij=x_i+\epsilon_2ij,\epsilon_1i|x_i~N(0,(2+0.2x_i)^2)}
#' \deqn{x_i\sim Uniform[10-40]}{x_i~Uniform[10-40]}
#'
#' for \eqn{i = 1, 2, \ldots, 100}, \eqn{j=1,2,\ldots,n_{2i}}{j=1,2,\ldots,n_2i} and
#'  the number of repeated measurements for each subject \eqn{i}
#' from the reference standard was \eqn{n_{2i} \sim Uniform[10-15]}{n_2i ~ Uniform[10-15]}.
#'
#'
"data1"


#' Simulated dataset 2
#'
#' In the simulated dataset 2, each subject has 1 to 5 measurement values from the new method
#' and 10 to 15 measurement values from the reference method.
#' Compared to the reference method, the new method has
#' differential bias of -4 and proportional bias of 1.2. Variance of the new method
#' is smaller than that for the referene method.
#'
#'  @format A data frame with three variables:
#' \describe{
#' \item{\code{id}}{identification number for subjects}
#' \item{\code{y1}}{values from the new measurement method}
#' \item{\code{y2}}{values from the reference measurement method}
#' }
#' @details
#' Dataset 2 was created based on the following equations:
#' \deqn{y_{1ij}=-4+1.2x_{i}+\varepsilon_{1ij}, \varepsilon_{1ij} \mid x_i \sim N(0,(1+0.1x_i)^2)}{
#' y_1ij=-4+1.2*x_i+\epsilon_1ij,\epsilon_1ij|x_i~N(0,1+0.1*x_i)^2)}
#' \deqn{y_{2ij}=x_i+\varepsilon_{2ij},\varepsilon_{2ij} \mid x_i\sim N(0,(2+0.2x_i)^2)}{
#' y_2ij=x_i+\epsilon_2ij,\epsilon_1i|x_i~N(0,(2+0.2x_i)^2)}
#' \deqn{x_i\sim Uniform[10-40]}{x_i~Uniform[10-40]}
#'
#' for \eqn{i = 1, 2, \ldots, 100}, \eqn{j=1,2,\ldots,n_{1i} / n_{2i}}{j=1,2,\ldots,n_1i or n_2i} and
#'  the number of repeated measurements for each subject \eqn{i}
#' from the new and reference method was \eqn{n_{1i} \sim Uniform[1-5]}{n_1i ~ Uniform[1-5]} and
#' \eqn{n_{2i} \sim Uniform[10-15]}{n_2i ~ Uniform[10-15]} respectively.
#'
#'
"data2"


#' Simulated dataset 3
#'
#' In the simulated dataset 3, each subject has 1 to 5 measurement values from the new method
#' and 10 to 15 measurement values from the reference method.
#' Compared to the reference method, the new method has
#' differential bias of 3 and proportional bias of 0.9. Variance of the new method
#' is larger than that for the referene method.
#'
#'  @format A data frame with three variables:
#' \describe{
#' \item{\code{id}}{identification number for subjects}
#' \item{\code{y1}}{values from the new measurement method}
#' \item{\code{y2}}{values from the reference measurement method}
#' }
#' @details
#' Dataset 3 was created based on the following equations:
#' \deqn{y_{1ij}=3+0.9x_{i}+\varepsilon_{1ij}, \varepsilon_{1ij} \mid x_i \sim N(0,(2+0.06x_i)^2)}{
#' y_1ij=3+0.9*x_i+\epsilon_1ij,\epsilon_1ij|x_i~N(0,2+0.06*x_i)^2)}
#' \deqn{y_{2ij}=x_i+\varepsilon_{2ij},\varepsilon_{2ij} \mid x_i\sim N(0,(1+0.01x_i)^2)}{
#' y_2ij=x_i+\epsilon_2ij,\epsilon_1i|x_i~N(0,(1+0.01x_i)^2)}
#' \deqn{x_i\sim Uniform[10-40]}{x_i~Uniform[10-40]}
#'
#' for \eqn{i = 1, 2, \ldots, 100}, \eqn{j=1,2,\ldots,n_{1i} / n_{2i}}{j=1,2,\ldots,n_1i or n_2i} and
#'  the number of repeated measurements for each subject \eqn{i}
#' from the new and reference method was \eqn{n_{1i} \sim Uniform[1-5]}{n_1i ~ Uniform[1-5]} and
#' \eqn{n_{2i} \sim Uniform[10-15]}{n_2i ~ Uniform[10-15]} respectively.
#'
#'
"data3"



