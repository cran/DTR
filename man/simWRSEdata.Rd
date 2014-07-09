\name{simWRSEdata}
\alias{simWRSEdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate data from sequentially randomized clinical trials (Guo and Tsiatis 2005)
}
\description{
This function generates a data set from sequentially randomized clinical trials as described in the simulation study of Guo and Tsiatis (2005). Because different assignments at the first stage are independent to each other,  the function only generates data for one of the assignments at the first stage. For example, if there are two first-stage therapies A1 and A2 available, the function only simulates the data for one of the two arms (e.g. A1). \cr

}
\usage{
simWRSEdata(n, max.c, pi.r, pi.z, mean.T0, 
mean.TR, mean.T1, mean.T2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
number of subjects assigned to A1
}
  \item{max.c}{
%%     ~~Describe \code{max.c} here~~
censoring time C is generated from uniform(0, max.c)
}
  \item{pi.r}{
%%     ~~Describe \code{pi.r} here~~
the probability of response. Response status R is generated from Bernoulli(pi.r)
}
  \item{pi.z}{
%%     ~~Describe \code{pi.z} here~~
the probability of being assigned to B2 among responders. The second-stage treatment indicator Z is generated from Bernoulli(pi.z) among responders (R=1). Z=0 if assigned to B1 at the second stage, and Z=1 if assigned to B2 at the second stage. For non-responders (R=0), Z=0
}
  \item{mean.T0}{
%%     ~~Describe \code{mean.T0} here~~
for nonresponders (R=0), a survival time T0 is drawn from exponential(1/mean.T0) with mean equal to mean.T0
}
  \item{mean.TR}{
%%     ~~Describe \code{mean.TR} here~~
for responders (R=1), a response time TR is drawn from exponential(1/mean.TR) with mean equal to mean.TR
}
  \item{mean.T1}{
%%     ~~Describe \code{mean.T1} here~~
the survival time from the response/consent to event time if assigned to B1, T1* is generated from exponential(1/mean.T1) with mean equal to mean.T1
}
  \item{mean.T2}{
%%     ~~Describe \code{mean.T2} here~~
the survival time from the response/consent to event time if assigned to B2, T2* is generated from exponential(1/mean.T2) with mean equal to mean.T2
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2.
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: TR, R, Z, U, and delta.
\item{TR}{a time to response for responders (R=1)}
\item{R}{response status, R=1 for responders, and R=0 for non-responders}
\item{Z}{second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2}
\item{U}{observed survival time, U is event time if delta=1, and U is censoring time if delta=0}
\item{delta}{censoring indicator, delta=1 for event, and delta=0 for censored}
}
\references{
%% ~put references to the literature/web site here ~
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
}
\note{
%%  ~~further notes~~
Set a seed right before \code{\link{simWRSEdata}} will help reproduce the same data.\cr
}
\seealso{
\code{\link{WRSEestimate}}
}
\examples{
 
#-------Example function usage-----------
n <- 100
max.c <- 3.5*365 
pi.r <- 0.5 
pi.z <- 0.5 
mean.T0 <- 182.5 
mean.TR <- 365 
mean.T1 <- 365 
mean.T2 <- 547.5

# Generate full data from SRD
set.seed(123)
data.A1 <- simWRSEdata(n,max.c,pi.r,pi.z,
mean.T0,mean.TR,mean.T1,mean.T2)
data.A2 <- simWRSEdata(n,max.c,pi.r,pi.z,
mean.T0,mean.TR,mean.T1,mean.T2)
WRSEdata <- cbind(X=c(rep(0,n), rep(1,n)), 
rbind(data.A1, data.A2)) # X=0 for A1; X=1 for A2
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
