\docType{data}
\name{data1}
\alias{data1}
\title{Simulated dataset 1}
\format{An object of class \code{data.frame} with 1255 rows and 3 columns.}
\usage{
data1
}
\description{
In the simulated dataset 1, each subject has 1 measurement value from the new method
and 10 to 15 measurement values from the reference method.
Compared to the reference method, the new method has
differential bias of -4 and proportional bias of 1.2. Variance of the new method
is smaller than that for the reference method.
}
\details{
A data frame with three variables:
\describe{
\item{\code{id}}{identification number for subjects}
\item{\code{y1}}{values from the new measurement method}
\item{\code{y2}}{values from the reference measurement method}
}

Dataset 1 was created based on the following equations:
\deqn{y_{1i}=-4+1.2x_i+\varepsilon_{1i}, \varepsilon_{1i} \mid x_i \sim N(0,(1+0.1x_i)^2)}{
y_1i=-4+1.2*x_i+\epsilon_1i,	\epsilon_1i | x_i ~ N(0,1+0.1*x_i)^2)}
\deqn{y_{2ij}=x_i+\varepsilon_{2ij},\varepsilon_{2ij} \mid x_i\sim N(0,(2+0.2x_i)^2)}{
y_2ij=x_i+\epsilon_2ij,	\epsilon_1i | x_i ~ N(0,(2+0.2x_i)^2)}
\deqn{x_i\sim Uniform[10-40]}{x_i ~ Uniform[10-40]}

for \eqn{i = 1, 2, \ldots, 100}, \eqn{j=1,2,\ldots,n_{2i}}{j=1,2,\ldots,n_2i} and
 the number of repeated measurements for each subject \eqn{i}
from the reference standard was \eqn{n_{2i} \sim Uniform[10-15]}{n_2i ~ Uniform[10-15]}.
}
\keyword{datasets}

