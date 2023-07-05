\name{YPmodelPhreg-package}
\alias{YPmodelPhreg-package}

\title{
The Short-Term and Long-Term Hazard Ratio Model with Proportional Adjustment
}
\description{
This package fits a semiparametric short term and long term hazard ratio model for the group variable, adjusted for proportional covariate effects. The model was studied in Yang and Prentice (2015), and it extends the two sample results in Yang and Prentice (2005). The model accommodates various non-proportional hazards scenarios including those with crossing hazards or survival functions. The main function of the package is \code{\link{ypreg}}, which provides the estimates and confidence intervals  of the model parameters, the average hazard ratio, and the hazard ratios at user-specified time points. Plots of point-wise and simultaneous confidence bands are also given.


}

\details{
\tabular{ll}{
Package: \tab YPmodelPhreg\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2023-04-21\cr
License: \tab GPL (>= 3)\cr
}
}

\author{
Daewoo Pak and Song Yang
Maintainer: Daewoo Pak <dpak@yonsei.ac.kr>
}

\value{
No return value, called for side effects
}

\references{
Yang, S., & Prentice, R. (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. \emph{Biometrika}, 92(1), 1-17.

Yang, S., & Prentice, R. L. (2015). Assessing potentially time-dependent treatment effect from clinical trials and observational studies for survival data, with applications to the Women's Health Initiative combined hormone therapy trial. \emph{Statistics in medicine}, 34(11), 1801-1817.
}
\keyword{ package }
\seealso{
  \code{\link{ypreg}}
}
\examples{
library(YPmodelPhreg)
data(colonexample)
head(colonexample)

res <- ypreg(colonexample, time.hr = c(1, 7))
res
plot(res)
}
