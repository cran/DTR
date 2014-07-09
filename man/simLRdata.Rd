\name{simLRdata}
\alias{simLRdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate data from sequentially randomized clinical trials (Kidwell and Wahed 2013)
}
\description{
This function generates a data set from sequentially randomized clinical trials as described in the simulation study of Kidwell and Wahed (2013). 
}
\usage{
simLRdata(n,max.c,pi.x,pi.r,pi.z,
mean.NR.1,mean.NR.2,mean.R.1,mean.R.2,
mean.RE.11,mean.RE.12,mean.RE.21,mean.RE.22)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
total number of subjects participating in the clinical trial
}
  \item{max.c}{
%%     ~~Describe \code{max.c} here~~
censoring time C is generated from uniform(0, max.c)
}
  \item{pi.x}{
%%     ~~Describe \code{pi.x} here~~
the probability of being assigned to A2 at the first stage. The first-stage treatment indicator X is generated from Bernoulli(pi.x). X=0 if assigned to A1 at the first stage, and X=1 if assigned to A2 at the first stage.
}
  \item{pi.r}{
%%     ~~Describe \code{pi.r} here~~
the probability of response. Response status R is generated from Bernoulli(pi.r)
}
  \item{pi.z}{
%%     ~~Describe \code{pi.z} here~~
the probability of being assigned to B2 among responders. The second-stage treatment indicator Z is generated from Bernoulli(pi.z) among responders (R=1). Z=0 if assigned to B1 at the second stage, and Z=1 if assigned to B2 at the second stage. For non-responders (R=0), Z=0
}
  \item{mean.NR.1}{
%%     ~~Describe \code{mean.NR.1} here~~
for patients who are assigned to A1 (X=0) and do not respond (R=0), a survival time T.NR.1 is drawn from exponential(1/mean.NR.1) with mean equal to mean.NR.1
}
  \item{mean.NR.2}{
%%     ~~Describe \code{mean.NR.2} here~~
for patients who are assigned to A2 (X=1) and do not respond (R=0), a survival time T.NR.2 is drawn from exponential(1/mean.NR.2) with mean equal to mean.NR.2
}
  \item{mean.R.1}{
%%     ~~Describe \code{mean.R.1} here~~
for patients who are assigned to A1 (X=0) and respond (R=1), time to response T.R.1 is drawn from exponential(1/mean.R.1) with mean equal to mean.R.1
}
  \item{mean.R.2}{
%%     ~~Describe \code{mean.R.2} here~~
for patients who are assigned to A2 (X=1) and respond (R=1), time to response T.R.2 is drawn from exponential(1/mean.R.2) with mean equal to mean.R.2
}
  \item{mean.RE.11}{
%%     ~~Describe \code{mean.RE.11} here~~
for patients who are assigned to A1 (X=0), respond (R=1), and then assigned to B1 (Z=0), a time from response to event is generated from exponential(1/mean.RE.11) with mean equal to mean.RE.11
}
  \item{mean.RE.12}{
%%     ~~Describe \code{mean.RE.12} here~~
for patients who are assigned to A1 (X=0), respond (R=1), and then assigned to B2 (Z=1), a time from response to event is generated from exponential(1/mean.RE.12) with mean equal to mean.RE.12
}
  \item{mean.RE.21}{
%%     ~~Describe \code{mean.RE.21} here~~
for patients who are assigned to A2 (X=1), respond (R=1), and then assigned to B1 (Z=0), a time from response to event is generated from exponential(1/mean.RE.21) with mean equal to mean.RE.21
}
  \item{mean.RE.22}{
%%     ~~Describe \code{mean.RE.22} here~~
for patients who are assigned to A2 (X=1), respond (R=1), and then assigned to B2 (Z=1), a time from response to event is generated from exponential(1/mean.RE.22) with mean equal to mean.RE.22
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2.
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: X, TR, R, Z, U, and delta.
\item{X}{first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2}
\item{TR}{time to response for responders (R=1)}
\item{R}{response status, R=1 for responders, and R=0 for non-responders}
\item{Z}{second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2}
\item{U}{observed survival time, U is event time if delta=1, and U is censoring time if delta=0}
\item{delta}{censoring indicator, delta=1 for event, and delta=0 for censored}
}
\references{
%% ~put references to the literature/web site here ~
Kidwell KM, Wahed AS: Weighted log-rank statistic to compare shared-path. adaptive treatment strategies. Biostatistics. 14(2):299-312, 2013
}
\note{
%%  ~~further notes~~
Set a seed right before \code{\link{simLRdata}} will help reproduce the same data.\cr
}
\seealso{
\code{\link{contrast_logrank}}
}
\examples{
 
#-------Example function usage-----------
n <- 100
max.c <- 12 # 30\% censoring
pi.x <- 0.5
pi.r <- 0.6 # 60\% response rate
pi.z <- 0.5
mean.NR.1 <- 1
mean.NR.2 <- 1
mean.R.1 <- 1
mean.R.2 <- 1
mean.RE.11 <- 5
mean.RE.12 <- 5
mean.RE.21 <- 5
mean.RE.22 <- 5

set.seed(123)
LRdata <- simLRdata(n,max.c,pi.x,pi.r,
pi.z,mean.NR.1,mean.NR.2,mean.R.1,mean.R.2,
mean.RE.11,mean.RE.12,mean.RE.21,mean.RE.22)
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
