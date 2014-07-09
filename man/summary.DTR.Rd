\name{summary.DTR}
\alias{summary.DTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of survival curves 
}
\description{
Returns an object of class \code{summary.DTR}. See \code{DTR.object} for details.
}
\usage{
\method{summary}{DTR}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{est} here~~
the result of a call to the \code{\link{LDTestimate}} function or \code{\link{WRSEestimate}} function
}
\item{\dots}{for future methods}
}
\value{
The function returns an object of class \code{summary.DTR}. \cr
}
\seealso{
\code{\link{DTR.object}}, \code{\link{print.DTR}},
\code{\link{print.summary.DTR}}, \code{\link{plot.DTR}} 
}
\examples{
\dontrun{
data("LDTdata")
est <- LDTestimate(data=LDTdata)
summary(est)}
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

