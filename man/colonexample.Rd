\name{colonexample}
\alias{colonexample}
\docType{data}
\title{An example from Chemotherapy for Stage B/C colon cancer}
\description{
The example \code{colonexample} in this package was created using data from one of the first successful trials of adjuvant chemotherapy for colon cancer. The original data set is available from \code{colon} in the \code{survival} library. The data set \code{colonexample} here uses only records with \code{etype = 2} and \code{rx = 1 or 3}, and it includes three variables: \code{extent}, \code{surg}, and \code{node4}.

}
\usage{data("colonexample")}
\format{
  A data frame with 619 observations on the following 6 variables.
  \describe{
    \item{\code{time}}{time until event or censoring (in years)}
    \item{\code{event}}{censoring status (1 = event, 0 = censored)}
    \item{\code{group}}{group indicator (1 = Lev(amisole)+5-FU), 0 = Obs(ervation))}
    \item{\code{X1}}{a binary covariate for \code{extent} (1 = serosa or contiguous structures, 0 = submucosa or muscle)}
    \item{\code{X2}}{a binary covariate for \code{surg} (1=long, 0=short)}
    \item{\code{X3}}{a binary covariate for \code{node4} (1 = more than 4 positive lymph nodes, 0 = otherwise)}
  }
}
\examples{
library(YPmodelPhreg)
data(colonexample)
head(colonexample)
}
\seealso{\code{\link{ypreg}}}
\keyword{datasets}
