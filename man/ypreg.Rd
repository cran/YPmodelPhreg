\name{ypreg}
\alias{ypreg}
\alias{ypreg.default}
\alias{print.ypreg}
\alias{plot.ypreg}
\title{
Fit a short-term and long-term hazard ratio model with proportional adjustment
}
\description{

The main results of the function are the estimations of:
    \itemize{
      \item parameters in the covariate-adjusted short-term and long-term hazard ratio model with confidence intervals;
      \item the average hazard ratio with the confidence interval; and
      \item the hazard ratio function along with point-wise and simultaneous confidence bands (confidence intervals for the hazard ratios at specific user input time points are also given).
    }
}

\usage{
\method{ypreg}{default}(data, alpha = 0.05, time.hr = NULL,
  L = NULL, U = NULL, repnum = 5000, tau = NULL, ...)
}

\arguments{
\item{...}{for S4 method only.}
\item{data}{A numeric matrix containing all variables in the data set. The columns must follow this order: 1) time until event or censoring, 2) censoring status (1 = event, 0 = censored), 3) binary group indicator taking values of 0 and 1 (e.g., 1 = treatment, 0 = control for a randomized trial), and 4) a set of numeric vectors of covariates. See the data structure of \code{\link{colonexample}}.}
\item{alpha}{A numeric value for the significance level. The default is \eqn{0.05}.}
\item{time.hr}{A numeric vector of time points at which hazard ratios will be estimated along with confidence intervals.}
\item{L}{A numeric value for the lower bound of the range [\code{L}, \code{U}] over which the simultaneous confidence bands for the hazard function are calculated. There must be at least one uncensored observation after the value of \code{L} and \code{L} < \code{U}.
}
\item{U}{A numeric value for the upper bound of the range [\code{L}, \code{U}] over which the simultaneous confidence bands for the hazard function are calculated. There must be at least one uncensored observation before the value of \code{U} and \code{U} > \code{L}.
}
\item{repnum}{The number of replications for the re-sampling method. The default is \eqn{5000}.}
\item{tau}{A numeric value for the maximum follow-up time. The default is \eqn{0.9\times}the maximum of all observations.}
}

\details{
The confidence intervals for the hazard ratios are obtained using the logarithmic transformation. When the user input interval [\code{L}, \code{U}] is different from the default interval, the intersection of the user input interval and the default interval is used. The point-wise confidence intervals and the simultaneous confidence bands can be plotted by supplying the object being returned by the function \code{\link{ypreg}} to the function \code{\link{plot.ypreg}}.
}


\value{
an object of S3 \code{ypreg} class representing the fit. The object also includes the results of the Cox proportional hazards model, implemented by using the \code{coxph} function in the \code{survival} library.

A list with at least the following elements:
\item{fit_coxph}{estimation results from the Cox proportional hazards model}
\item{best_b0}{the estimates from the short-term and long-term hazard ratio model without proportional adjustment}
\item{best_ypx}{the estimates from the short-term and long-term hazard ratio model with proportional adjustment}
\item{res_summ}{summary of estimation results with the covariate-adjusted short-term and long-term hazard ratio model}
\item{res_hrci}{estimation results of hazard ratios at \code{time.hr}}
}

\references{
Yang, S., & Prentice, R. (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. \emph{Biometrika}, 92(1), 1-17.

Yang, S., & Prentice, R. L. (2015). Assessing potentially time-dependent treatment effect from clinical trials and observational studies for survival data, with applications to the Women's Health Initiative combined hormone therapy trial. \emph{Statistics in medicine}, 34(11), 1801-1817.
}

\keyword{ypreg}
\seealso{
  \code{\link{plot.ypreg}}
}

\examples{
library(YPmodelPhreg)
data(colonexample)
head(colonexample)

res <- ypreg(colonexample, time.hr = c(1, 7))
res
plot(res)

}


