\name{bland_altman_plot}
\alias{bland_altman_plot}
\title{Extended Bland-Atlman limits of agreement (LoA) plot}
\usage{
bland_altman_plot(data, new = "y1", Ref = "y2", ID = "id", fill = TRUE)
}
\arguments{
\item{data}{a dataframe contains the object identification number (id),
the measurement values from the new measurement method (y1) and those from
the reference standard (y2)}

\item{new}{specify the variable name or location for the new measurement method}

\item{Ref}{specify the variable name or location for the reference  measuerment method}

\item{ID}{specify the variable name for location for the subject identification number (id)}

\item{fill}{logical. if \code{TRUE} use the average value for new methods to
fill out the missing value (only useful for drawing a plot with all the
measurements by the reference standard)}
}
\description{
This function produces the extended Bland-Altman Limits of Agreement (LoA) plot
when there are repeated measurements with possibly heteroscedastic variance of 
measurement errors.
}
\details{
This functions computes the limits of agreement (LoA) when there are
repeated measurements and possibly heteroscedastic variance of measurement errors
}
\examples{
### Load the data
data(data1)
### Bland and Altman's plot
bland_altman_plot(data1)
}
\author{
Mingkai Peng & Patrick Taffé
}

