\name{print.DTR}
\alias{print.DTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a short summary of survival curves
}
\description{
Print number of observations, number of events, and the median survival with confidence limits for the  
median for dynamic treatment regimes 
}
\usage{
\method{print}{DTR}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"DTR"}, which is the result of the \code{LDTestimate} function or \code{WRSEestimate} function }
\item{\dots}{for future results}
}
\details{
The median and its confidence interval are defined by drawing a horizontal line at 0.5 on the plot of the survival curve and its confidence bands.
}
\seealso{
\code{\link{DTR.object}}, \code{\link{summary.DTR}}, 
\code{\link{print.summary.DTR}}, \code{\link{plot.DTR}} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
