\name{MethodCompare}
\alias{MethodCompare}
\title{MethodCompare}
\description{
The package "MethodCompare" allows one to assess bias and precision of a new measurement method with respect to a reference method (also called "reference standard").
It requires repeated measurements by at least one of the two measurement methods.  

In this implementation, it is assumed by default that the reference method has repeated measurements
and the new method may have as few as only one measurement per individual (The methodology can be adapted if you have more repeated measurements 
by the new method than by the reference method, see ref. below).  

A manuscript with details concerning the methodology and its application can be found:

<doi:10.1177/0962280218759693>

It implements the methodology developped in:  

Taffé P. Effective plots to assess bias and precision in method comparison studies. Stat Methods Med Res 2018;27:1650-1660.  

NB: Further methodological developpements have been made and will be implemented in a future version of the package:

Taffé P. Assessing bias, precision, and agreement in method comparison studies. Stat Methods Med Res 2020;29:778-796. doi:10.1177/0962280218759693

For other relevant references:

Taffé P, Peng M, Stagg V, Williamson T. Biasplot: A package to effective plots to assess bias and precision in method comparison studies. Stata J 2017;17:208-221.

Taffé P, Peng M, Stagg V, Williamson T. MethodCompare: An R package to assess bias and precision in method comparison studies. Stat Methods Med Res 2019;28:2557-2565.

Taffé P, Halfon P, Halfon M. A new statistical methodology to assess bias and precision overcomes the defects of the Bland & Altman method. J Clin Epidemiol 2020;124:1-7.

Taffé P. When can the Bland-Altman limits of agreement method be used and when it should not be used. J Clin Epidemiol 2021; 137:176-181.
}
\author{
Mingkai Peng & Patrick Taffé
}