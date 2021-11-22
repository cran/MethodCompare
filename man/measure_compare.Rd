\name{measure_compare}
\alias{measure_compare}
\title{Estimation of the amount of bias of the new measurement method relative to the
reference method, with possibly heteroscedastic variance of the measurement errors}
\usage{
measure_compare(data, new = "y1", Ref = "y2", ID = "id")
}
\arguments{
\item{data}{a dataframe containing the identification number of the subject (id),
the measurement values from the new measurement method (y1) and those
from the reference method).}

\item{new}{specify the variable name or location of the new measurement method}

\item{Ref}{specify the variable name or location of the reference standard}

\item{ID}{specify the variable name for location of the subject identification
number}
}
\value{
The function returns a list with the following items:

\itemize{
  \item Bias: differential and proportional bias for new method and the
  associated 95 percent confidence intervals
  \item Models: list of models fitted in estimation procedure
  \item Ref: a data frame containing the various variables used to
  compute the bias and precision plots, as well the smooth standard errors
  estimates of the reference standard
  \item New: a data frame containing the various variables used to compute
  the bias and precision plots, as well the smooth standard errors estimates
   of the new measurement method
}
}
\description{
This function implements the methodology reported in the paper:
Taffé P. Effective plots to assess bias and precision in method comparison studies. Stat Methods Med Res 2018;27:1650-1660.
Other relevant references:
Taffé P, Peng M, Stagg V, Williamson T. Biasplot: A package to effective plots to assess bias and precision in method comparison studies. Stata J 2017;17:208-221.
Taffé P, Peng M, Stagg V, Williamson T. MethodCompare: An R package to assess bias and precision in method comparison studies. Stat Methods Med Res 2019;28:2557-2565.
Taffé P, Halfon P, Halfon M. A new statistical methodology to assess bias and precision overcomes the defects of the Bland & Altman method. J Clin Epidemiol 2020;124:1-7.
Taffé P. Assessing bias, precision, and agreement in method comparison studies. Stat Methods Med Res 2020;29:778-796.
Taffé P. When can the Bland-Altman limits of agreement method be used and when it should not be used. J Clin Epidemiol 2021; 137:176-181.
}
\details{
This functions implements the new estimation procedure to assess bias and precision of a new measurement method
with respect to a reference standard, as well as
Bland & Altman's limits of agreement extended to the setting of
possibly heteroscedastic variance of the measurement errors.
}
\examples{
### Load the data
data(data1)
### Analysis
measure_model <- measure_compare(data1)
}
\author{
Mingkai Peng & Patrick Taffé
}

