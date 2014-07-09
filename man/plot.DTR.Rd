\name{plot.DTR}
\alias{plot.DTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for \code{DTR} objects 
}
\description{
Plot the survival estimates and their 95\% confidence bands for each dynamic treatment regime
}
\usage{
\method{plot}{DTR}(x, confidence.interval = FALSE, 
xlab = "Time", ylab = "Survival probability", 
line.color = c("black", "grey40", "grey60", "grey80"), 
legend.position = "right", censored=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
an object of class \code{DTR}, usually returned by the \code{LDTestimate} function or \code{WRSEestimate} function.
}
\item{confidence.interval}{
If confidence.interval=FALSE (default), the 95\% confidence bands are not plotted. If confidence.interval=TRUE, the 95\% confidence bands are plotted as shadows.
}
\item{xlab}{
label given to the x-axis. Default is "Time".
}
\item{ylab}{
label given to the y-axis. Default is "Survival probability".
}
\item{line.color}{
colors for the lines. Default are "black", "grey40", "grey60", and "grey80" for A1B1, A1B2, A2B1, and A2B2 respectively.
}
\item{legend.position}{
the position of legend: "left", "right" (default), "bottom", "top", or two-element numeric vector (e.g. c(0.6,0.9)) 
}
\item{censored}{
If censored=FALSE (default), the censoring ticks are not plotted. If censored=TRUE, the censoring times are plotted as ticks  
}
\item{\dots}{for future methods}
}
\seealso{
\code{\link{DTR.object}}, \code{\link{print.DTR}}, \code{\link{summary.DTR}}, \code{\link{print.summary.DTR}}
}
\examples{
\dontrun{
data("LDTdata")
est <- LDTestimate(data=LDTdata)
plot(est, confidence.interval=TRUE, censored=TRUE)

data("WRSEdata")
est <- WRSEestimate(data=WRSEdata)
plot(est)}
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
