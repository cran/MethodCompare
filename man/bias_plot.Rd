% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bias_plot.R
\name{bias_plot}
\alias{bias_plot}
\title{Plot the bias and measurements}
\usage{
bias_plot(object)
}
\arguments{
\item{object}{list returned by \link{measure_compare} function.}
}
\description{
This function draws the "bias plot", which is used to visually assess the
bias of the new method relative to the reference method. It is obtained by
graphing a scatter plot of \code{y1} (new method) and \code{y2} (reference method) versus
the BLUP of the latent trait, \code{x}, along with the two regression lines.
The function adds a second scale on the right axis, showing the relationship
between the estimated amount of bias and BLUP of the latent trait, \code{x}.
}
\examples{
\donttest{
### Load the data
data(data1)
### Analysis
measure_model <- measure_compare(data1, nb_simul=100)
### Plot the bias
bias_plot(measure_model)}
}
