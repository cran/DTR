\name{summary.CHR}
\alias{summary.CHR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of cumulative hazard ratio estimates 
}
\description{
Returns an object of class \code{summary.CHR}. See \code{CHR.object} for details.
}
\usage{
\method{summary}{CHR}(object, log.CHR=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{est} here~~
the result of a call to the \code{\link{CHRestimate}} function
}
  \item{log.CHR}{
%%     ~~Describe \code{est} here~~
if log.CHR=FALSE (default), the summary of cumulative hazard ratio estimates is returned; if log.CHR=TRUE, the summary of log cumulative hazard ratio estimates is returned
}
\item{\dots}{for future methods}
}
\value{
The function returns an object of class \code{summary.CHR}. \cr
}
\seealso{
\code{\link{CHR.object}}, \code{\link{print.CHR}},
\code{\link{print.summary.CHR}}, \code{\link{plot.CHR}} 
}
\examples{
\dontrun{
data("CHRdata")
est <- CHRestimate(data=CHRdata, covar="V1")
summary(est, log.CHR=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}

