\name{compare_plot}
\alias{compare_plot}
\title{Comparison plot}
\usage{
compare_plot(object)
}
\arguments{
\item{object}{an object returned by a call to \link{measure_compare}}
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
Mingkai Peng & Patrick Taffé
}

