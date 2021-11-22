\name{precision_plot}
\alias{precision_plot}
\title{Precision plot}
\usage{
precision_plot(object)
}
\arguments{
\item{object}{an object returned by a call to \link{measure_compare}}
}
\description{
This plot allows the visual comparison of the precision (i.e. standard deviation)
of the new measurement method with that of the reference standard, by creating
a scatter plot of the estimated standard deviations against the best linear
prediction (BLUP) of the true latent trait x.
}
\examples{
### load the data
data(data1)
### analysis
measure_model <- measure_compare(data1)
### Precision plot
precision_plot(measure_model)

}
\author{
Mingkai Peng & Patrick Taffé
}

