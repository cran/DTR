\name{plot.CHR}
\alias{plot.CHR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for \code{CHR} objects 
}
\description{
Plot the cumulative hazard ratio estimates and their 95\% confidence bands for the comparison between two dynamic treatment regimes
}
\usage{
\method{plot}{CHR}(x, log.CHR = FALSE, confidence.interval = FALSE, 
xlab = "Time", line.color = c("black", "grey30", "grey50", 
"grey60", "grey70", "grey80"), legend.position = "right", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
an object of class \code{CHR}, usually returned by the \code{CHRestimate} function.
}
\item{log.CHR}{
%%     ~~Describe \code{est} here~~
if log.CHR=FALSE (default), the cumulative hazard ratio estimates are plotted; if log.CHR=TRUE, the log cumulative hazard ratio estimates are plotted
}
\item{confidence.interval}{
If confidence.interval=FALSE (default), the 95\% confidence bands are not plotted. If confidence.interval=TRUE, the 95\% confidence bands are plotted as shadows.
}
\item{xlab}{
label given to the x-axis. Default is "Time".
}
\item{line.color}{
colors for the lines. Default are "black", "grey30", "grey50", "grey60", "grey70", and "grey80" for A1B2 vs. A1B1, A2B1 vs. A1B1, A2B2 vs. A1B1, A2B1 vs. A1B2, A2B2 vs. A1B2, and A2B2 vs. A2B1 respectively.
}
\item{legend.position}{
the position of legend: "left", "right" (default), "bottom", "top", or two-element numeric vector (e.g. c(0.6,0.9)) 
}
\item{\dots}{for future methods}
}
\seealso{
\code{\link{CHR.object}}, \code{\link{print.CHR}}, \code{\link{summary.CHR}}, \code{\link{print.summary.CHR}}
}
\examples{
\dontrun{
data("CHRdata")
est <- CHRestimate(data=CHRdata)
plot(est, confidence.interval=TRUE)
plot(est, log.CHR=TRUE, confidence.interval=FALSE)}
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
