\name{DTR.CHR.plot}
\alias{DTR.CHR.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cumulative hazard ratio/Log cumulative hazard ratio plot for dynamic treatment regimes (DTRs)
%%  ~~function to do ... ~~
}
\description{
This function plots the estimated cumulative hazard ratios (CHR) or the natural logarithms of the estimated CHRs between two different dynamic treatment regimes (DTRs) from sequentially randomized clinical trials using the CHR estimates proposed in Tang and Wahed (2013) [Epub ahead of print]. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is implemented in the current version of the function. \cr
}
\usage{
DTR.CHR.plot(fdata, log.CHR=FALSE, confidence.interval=FALSE, xlab="Time", 
line.color=c("black", "grey30", "grey50", "grey60", "grey70", "grey80"), 
legend.position="right")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdata}{
%%     ~~Describe \code{fdata} here~~
a data frame (X, R, Z, U, delta, V) representing the full data from a simple two stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored \cr
V: covariates to be adjusted. The function allows adjustment for one covariate or more than one covariates \cr
}
  \item{log.CHR}{
%%     ~~Describe \code{method} here~~
If log.CHR=FALSE (default), the estimated CHRs are plotted. If log.CHR=TRUE, the natural logarithms of estimated CHRs are plotted.
}
  \item{confidence.interval}{ 
%%     ~~Describe \code{confidence.interval} here~~  
If confidence.interval=FALSE (default), the 95\% confidence bands are not plotted. If confidence.interval=TRUE, the 95\% confidence bands are plotted as grey shadows.
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
label for x axis. Default is "Time".
}
  \item{line.color}{
%%     ~~Describe \code{line.color} here~~
colors for the lines. Default are "black", "grey30", "grey50", "grey60", "grey70" and "grey80" for A1B2 vs. A1B1, A2B1 vs. A1B1, A2B2 vs. A1B1, A2B1 vs. A1B2, A2B2 vs. A1B2, and A2B2 vs. A2B1 respectively.
}
  \item{legend.position}{
%%     ~~Describe \code{legend.position} here~~
the position of legend: "left", "right" (default), "bottom", "top", or two-element numeric vector (e.g. c(0.6,0.9)) 
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2. Based on four DTRs, six different cumulative hazard ratios (CHRs) are computed: CHR(A1B2 vs. A1B1), CHR(A2B1 vs. A1B1), CHR(A2B2 vs. A1B1), CHR(A2B1 vs. A1B2), CHR(A2B2 vs. A2B1), and CHR(A2B2 vs. A2B1). The natural logarithms of the CHRs are also computed.\cr
}
\references{
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print]
}
\seealso{
\code{\link{DTR}}, \code{\link{sim.CHR.data}}, \code{\link{updateBeta}}, \code{\link{CHR.estimator}}, \code{\link{CHR.Wald.test}}
}
\examples{
require(DTR)
  
#-------Example function usage-----------
n <- 800
max.c <- 5
pi.x <- 0.5
pi.r <- 0.6
pi.z <- 0.5
gamma10 <- 1
gamma11 <- 1.5
gamma12 <- 1.2
gamma20 <- 1
gamma21 <- 1.5
gamma22 <- 1.3
alpha10 <- 0.5
alpha11 <- 0.2
alpha12 <- 0.1
alpha20 <- 0.5
alpha21 <- 0.2
alpha22 <- 0.05
beta <- c(0.5, 0.5)

set.seed(123)
fdata <- sim.CHR.data(n,max.c,pi.x,pi.r,pi.z,gamma10,gamma11,gamma12,
gamma20,gamma21,gamma22,alpha10,alpha11,alpha12,alpha20,alpha21,alpha22,beta)

\dontrun{

# Plot CHR
DTR.CHR.plot(fdata, confidence.interval=TRUE)

# Plot log CHR
DTR.CHR.plot(fdata, log.CHR=TRUE, confidence.interval=FALSE)
}
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
