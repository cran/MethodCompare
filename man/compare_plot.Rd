% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Compare_plot.R
\name{compare_plot}
\alias{compare_plot}
\title{Plot used to visualize the recalibration of the new method after estimating the bias}
\usage{
compare_plot(object)
}
\arguments{
\item{object}{an object retunred by a call to \link{measure_compare}}
}
\description{
This function allows the visualization of the bias-corrected values (i.e.
recalibrated values, variable y1_corr) of the new measurement method.
}
\examples{
### load the data
data(data1)
### analysis
measure_model <- measure_compare(data1)
### compare plot
compare_plot(measure_model)

}
\author{
Mingkai Peng
}

