\name{DTR.surv.plot}
\alias{DTR.surv.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Survival plot for dynamic treatment regimes (DTRs)
%%  ~~function to do ... ~~
}
\description{
This function plots the survival functions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials based on the LDT estimator proposed in Lunceford, Davidian and Tsiatis (2002) or the WRSE estimator proposed in Guo and Tsiatis (2005). For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is implemented in the current version of the function. \cr
}
\usage{
DTR.surv.plot(fdata, method="LDT", L=.Machine$double.xmax, 
confidence.interval=FALSE, xlab="Time", ylab="Survival probability", 
line.color=c("black", "grey40", "grey60", "grey80"), 
legend.position="right")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdata}{
%%     ~~Describe \code{fdata} here~~
a data frame (X, TR, R, Z, U, delta) representing the full data from a simple two stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
TR: optional time to response for responders (R=1), only needed if method="WRSE" \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored 
}
  \item{method}{
%%     ~~Describe \code{method} here~~
the survival estimates to be used for the plot. The survival estimator proposed in Lunceford, Davidian and Tsiatis (2002) is used when method="LDT" (default), and the weighted risk set estimator (WRSE) proposed in Guo and Tsiatis (2005) is used when method="WRSE".
}
  \item{L}{
%%     ~~Describe \code{L} here~~
restricted survival time. Default is .Machine$double.xmax, which is the largest double value of R. Set L to a numeric number smaller than the maximum follow-up time if restricted follow-up time up to L is considered. The value of L is only used when method="LDT".
}
  \item{confidence.interval}{ 
%%     ~~Describe \code{confidence.interval} here~~  
If confidence.interval=FALSE (default), the 95\% confidence bands are not plotted. If confidence.interval=TRUE, the 95\% confidence bands are plotted as grey shadows.
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
label for x axis. Default is "Time".
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
label for y axis. Default is "Survival probability".
}
  \item{line.color}{
%%     ~~Describe \code{line.color} here~~
colors for the lines. Default are "black", "grey40", "grey60", and "grey80" for A1B1, A1B2, A2B1, and A2B2 respectively.
}
  \item{legend.position}{
%%     ~~Describe \code{legend.position} here~~
the position of legend: "left", "right" (default), "bottom", "top", or two-element numeric vector (e.g. c(0.6,0.9)) 
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2. \cr
}
\references{
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002 \cr
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
}
\seealso{
\code{\link{DTR}}, \code{\link{sim.LDT.data}}, \code{\link{LDT.estimator}}, \code{\link{sim.WRSE.data}}, \code{\link{WRSE.estimator}}
}
\examples{
require(DTR)
  
#-------Example function usage-----------
n <- 500
max.c <- 3.5*365 
pi.r <- 0.5 
pi.z <- 0.5 
mean.T0 <- 182.5 
mean.TR <- 365 
mean.T1 <- 365 
mean.T2 <- 547.5

# Generate A1 data from SRD
set.seed(123)  # Set a seed to reproduce the same data every time running the following command
data <- sim.WRSE.data(n,max.c,pi.r,pi.z,mean.T0,mean.TR,mean.T1,mean.T2)

# Generate full data from SRD
set.seed(123)
data.A1 <- sim.WRSE.data(n,max.c,pi.r,pi.z,mean.T0,mean.TR,mean.T1,mean.T2)
data.A2 <- sim.WRSE.data(n,max.c,pi.r,pi.z,mean.T0,mean.TR,mean.T1,mean.T2)
fdata <- cbind(X=c(rep(0,n), rep(1,n)), rbind(data.A1, data.A2)) # X=0 for A1; X=1 for A2

\dontrun{

# Plot survival function
DTR.surv.plot(fdata, method="LDT", confidence.interval=TRUE)
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
