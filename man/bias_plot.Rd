\name{bias_plot}
\alias{bias_plot}
\title{Bias plot}
\description{
This function draws the "bias plot", which is used to visually assess the bias of the new method relative to the reference method. 
It is obtained by graphing a scatter plot of y1 (new method) and y2 (reference method) versus the BLUP of y2 along with the two regression lines, and adds a second scale on the right axis showing the relationship between the estimated amount of bias and BLUP of y2.
}
\usage{
bias_plot(object)
}
\arguments{
\item{object}{an object returned by a call to \link{measure_compare}}
}
\examples{
### load the data
data(data1)
### analysis
measure_model <- measure_compare(data1)
### Bias plot
bias_plot(measure_model)
}
\author{
Mingkai Peng & Patrick Taffé
}