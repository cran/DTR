\name{print.CHR}
\alias{print.CHR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a short summary of cumulative hazard ratio estimates
}
\description{
Print the comparison between two dynamic treatment regimes, 75th percentile of the observed times, the cumulative hazard ratio estimates with 95\% confidence intervals at 75th percentile of the observed times for each comparison, and the log cumulative hazard ratio estimates with 95\% confidence intervals at 75th percentile of the observed times for each comparison.
}
\usage{
\method{print}{CHR}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"CHR"}, which is the result of the \code{CHRestimate} function}
\item{\dots}{for future results}
}
\details{
The 75th percentile of the observed times is selected.
}
\seealso{
\code{\link{CHR.object}}, \code{\link{summary.CHR}}, 
\code{\link{print.summary.CHR}}, \code{\link{plot.CHR}} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
