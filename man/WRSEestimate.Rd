\name{WRSEestimate}
\alias{WRSEestimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function for calculating WRSE estimates
}
\description{
The function computes the weighted risk set estimates (WRSE) of the survival functions and their estimated standard errors for dynamic treatment regimes (DTRs) at observed event times as proposed in Guo and Tsiatis (2005) Equation (3) and Equation (16).
}
\usage{
WRSEestimate(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame (X, TR, R, Z, U, delta) representing the data from a two-stage randomized designs with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: treatment arm, X=0 for A1, and X=1 for A2\cr
TR: the time to response for responders (R=1) \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0  \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2.
}
\value{
The function returns an object of class \code{DTR}. See \code{DTR.object} for details. \cr
}
\references{
%% ~put references to the literature/web site here ~
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
}
\seealso{
\code{\link{simWRSEdata}}, \code{\link{DTR.object}}, \code{\link{print.DTR}}, 
\code{\link{summary.DTR}}, \code{\link{print.summary.DTR}}, \code{\link{plot.DTR}}
}
\examples{
\dontrun{
data("WRSEdata")
est <- WRSEestimate(data=WRSEdata)
est}
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


