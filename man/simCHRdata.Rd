\name{simCHRdata}
\alias{simCHRdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate data from sequentially randomized clinical trials (Tang and Wahed 2013 [Epub ahead of print])
}
\description{
This function generates a data set from sequentially randomized clinical trials as described in the simulation study of Tang and Wahed (2013) [Epub ahead of print]. 
}
\usage{
simCHRdata(n,max.c,pi.x,pi.r,pi.z,gamma10,
gamma11,gamma12,gamma20,gamma21,gamma22,alpha10,
alpha11,alpha12,alpha20,alpha21,alpha22,beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
total number of subjects participating in the clinical trial
}
  \item{max.c}{
%%     ~~Describe \code{max.c} here~~
censoring time C is generated from uniform(max.c/2, max.c)
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
  \item{gamma10}{
%%     ~~Describe \code{gamma10 alpha10} here~~
for individuals who are assigned to A1 at first-stage, and do not respond, survival times are drawn from Weibull distribution with parameters alpha10 and gamma10
}
  \item{gamma11}{
%%     ~~Describe \code{gamma11 alpha11} here~~
for individuals who are assigned to A1 at first-stage, respond to A1, and are assigned to B1 at second-stage, survival times are drawn from Weibull distribution with parameters alpha11 and gamma11
}
  \item{gamma12}{
%%     ~~Describe \code{gamma12 alpha12} here~~
for individuals who are assigned to A1 at first-stage, respond to A1, and are assigned to B2 at second-stage, survival times are drawn from Weibull distribution with parameters alpha12 and gamma12
}
  \item{gamma20}{
%%     ~~Describe \code{gamma20 alpha20} here~~
for individuals who are assigned to A2 at first-stage, and do not respond, survival times are drawn from Weibull distribution with parameters alpha20 and gamma20
}
  \item{gamma21}{
%%     ~~Describe \code{gamma21 alpha21} here~~
for individuals who are assigned to A2 at first-stage, respond to A2, and are assigned to B1 at second-stage, survival times are drawn from Weibull distribution with parameters alpha21 and gamma21
}
  \item{gamma22}{
%%     ~~Describe \code{gamma22 alpha22} here~~
for individuals who are assigned to A2 at first-stage, respond to A2, and are assigned to B2 at second-stage, survival times are drawn from Weibull distribution with parameters alpha22 and gamma22
}
  \item{alpha10}{
%%     ~~Describe \code{gamma10 alpha10} here~~
for individuals who are assigned to A1 at first-stage, and do not respond, survival times are drawn from Weibull distribution with parameters alpha10 and gamma10
}
  \item{alpha11}{
%%     ~~Describe \code{gamma11 alpha11} here~~
for individuals who are assigned to A1 at first-stage, respond to A1, and are assigned to B1 at second-stage, survival times are drawn from Weibull distribution with parameters alpha11 and gamma11
}
  \item{alpha12}{
%%     ~~Describe \code{gamma12 alpha12} here~~
for individuals who are assigned to A1 at first-stage, respond to A1, and are assigned to B2 at second-stage, survival times are drawn from Weibull distribution with parameters alpha12 and gamma12
}
  \item{alpha20}{
%%     ~~Describe \code{gamma20 alpha20} here~~
for individuals who are assigned to A2 at first-stage, and do not respond, survival times are drawn from Weibull distribution with parameters alpha20 and gamma20
}
  \item{alpha21}{
%%     ~~Describe \code{gamma21 alpha21} here~~
for individuals who are assigned to A2 at first-stage, respond to A2, and are assigned to B1 at second-stage, survival times are drawn from Weibull distribution with parameters alpha21 and gamma21
}
  \item{alpha22}{
%%     ~~Describe \code{gamma22 alpha22} here~~
for individuals who are assigned to A2 at first-stage, respond to A2, and are assigned to B2 at second-stage, survival times are drawn from Weibull distribution with parameters alpha22 and gamma22
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
the coefficient vector for two covariates V1 and V2, for example, beta=c(0.5,0.5). Both covariates V1 and V2 are generated from Bernoulli(0.5)
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2.
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: X, R, Z, U, delta, V1 and V2.
\item{X}{first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2}
\item{R}{response status, R=1 for responders, and R=0 for non-responders}
\item{Z}{second-stage indicator among responders (R=1), Z=0 if assigned to B1, and Z=1 if assigned to B2}
\item{U}{observed survival time, U is event time if delta=1, and U is censoring time if delta=0}
\item{delta}{censoring indicator, delta=1 for event, and delta=0 for censored}
\item{V1}{a binary covariate following Bernoulli(0.5)}
\item{V2}{a binary covariate following Bernoulli(0.5)}
}
\references{
%% ~put references to the literature/web site here ~
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print]
}
\note{
%%  ~~further notes~~
Set a seed right before \code{\link{simCHRdata}} will help reproduce the same data.\cr
}
\seealso{
\code{\link{CHRestimate}}
}
\examples{

#-------Example function usage-----------
n <- 200
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
CHRdata <- simCHRdata(n,max.c,pi.x,pi.r,
pi.z,gamma10,gamma11,gamma12,gamma20,gamma21,
gamma22,alpha10,alpha11,alpha12,alpha20,alpha21,
alpha22,beta)
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
