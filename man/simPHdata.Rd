\name{simPHdata}
\alias{simPHdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate data from sequentially randomized clinical trials (Tang and Wahed 2011)
}
\description{
This function generates a data set from sequentially randomized clinical trials as described in the simulation study of Tang and Wahed (2011). 
}
\usage{
simPHdata(n, max.c, pi.x, pi.z, lambda, alpha, 
beta1, beta2, beta3, beta4, beta5, gamma)
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
  \item{pi.z}{
%%     ~~Describe \code{pi.z} here~~
the probability of being assigned to B2 among responders. The second-stage treatment indicator Z is generated from Bernoulli(pi.z) among responders (R=1). Z=0 if assigned to B1 at the second stage, and Z=1 if assigned to B2 at the second stage. For non-responders (R=0), Z=0
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
baseline hazard
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
a response time TR is drawn from exponential(alpha) with mean equal to 1/alpha
}
  \item{beta1}{
%%     ~~Describe \code{beta1} here~~
coefficient for first-stage indicator X
}
  \item{beta2}{
%%     ~~Describe \code{beta2} here~~
coefficient for time-varying response indicator R(t)=R*I(TR<t)
}
  \item{beta3}{
%%     ~~Describe \code{beta3} here~~
coefficient for the interaction between X and R(t)
}
  \item{beta4}{
%%     ~~Describe \code{beta4} here~~
coefficient for the interaction between R(t) and second-stage indicator Z
}
  \item{beta5}{
%%     ~~Describe \code{beta5} here~~
coefficient for the three-way interaction among X, R(t), and Z
}
  \item{gamma}{
%%     ~~Describe \code{beta2} here~~
coefficient for the covariate V. Covariate V is drawn from normal distribution with mean 1 and standard deviation 0.5
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2.
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: X, TR, R, Z, U, delta, and V.
\item{X}{first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2}
\item{TR}{time to response}
\item{R}{response status, R=1 for responders, and R=0 for non-responders}
\item{Z}{second-stage indicator among responders (R=1), Z=0 if assigned to B1, and Z=1 if assigned to B2}
\item{U}{observed survival time, U is event time if delta=1, and U is censoring time if delta=0}
\item{delta}{censoring indicator, delta=1 for event, and delta=0 for censored}
\item{V}{a continuous covariate following normal distribution with mean 1 and standard deviation 0.5}
}
\references{
%% ~put references to the literature/web site here ~
Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary variables. Journal of Applied Statistics 38(12):2925-2938, 2011
}
\note{
%%  ~~further notes~~
Set a seed right before \code{\link{simPHdata}} will help reproduce the same data.\cr
}
\seealso{
\code{\link{PHfit}}, \code{\link{contrast_ph}}
}
\examples{

#-------Example function usage-----------
n <- 400
pi.x <- 0.5
pi.z <- 0.5
lambda <- 1/4
alpha <- 1/6
beta1 <- -0.5
beta2 <- -0.8
beta3 <- 0.5
beta4 <- 1
beta5 <- -1
gamma <- -0.5
max.c <- 14

set.seed(123)
PHdata <- simPHdata(n,max.c,pi.x,pi.z,lambda,alpha,
beta1,beta2,beta3,beta4,beta5,gamma)
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
