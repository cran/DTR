\name{simLDTdata}
\alias{simLDTdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to simulate data from sequentially randomized clinical trials (Lunceford, Davidian and Tsiatis 2002)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates a data set from sequentially randomized clinical trials as described in the simulation study of Lunceford, Davidian and Tsiatis (2002). Because different assignments at the first stage are independent to each other,  the function only generates data for one of the assignments at the first stage. For example, if there are two first-stage therapies A1 and A2 available, the function only simulates the data for one of the two arms (e.g. A1). \cr

}
\usage{
simLDTdata(n, max.c, pi.r, pi.z, lambda, alpha, 
beta1, beta2, L=.Machine$double.xmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
number of subjects assigned to A1
}
  \item{max.c}{
%%     ~~Describe \code{c} here~~
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
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
for nonresponders (R=0), a survival time T*_lambda is drawn from exponential(lambda) with mean 1/lambda
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
for responders, a response time T*_alpha is drawn from exponential(alpha) with mean 1/alpha
}
  \item{beta1}{
%%     ~~Describe \code{beta1} here~~
post-response survival time under B1, T*_11 is drawn from exponential(e^beta1)
}
  \item{beta2}{
%%     ~~Describe \code{beta2} here~~
post-response survival time under B2, T*_12 is drawn from exponential(e^(beta1+beta2*T*_11))
}
  \item{L}{restricted survival time. Default is .Machine$double.xmax, which is the largest double value of R. Set L to a numeric number smaller than the maximum follow-up time if restricted follow-up time up to L is considered.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2.
}
\value{
The function returns a data set with columns: R, Z, U, and delta.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{R}{response status, R=1 for responders, and R=0 for non-responders}
\item{Z}{second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2}
\item{U}{observed survival time, U is event time if delta=1, and U is censoring time if delta=0}
\item{delta}{censoring indicator, delta=1 for event, and delta=0 for censored}
}
\references{
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
}
\note{
%%  ~~further notes~~
Set a seed right before \code{\link{simLDTdata}} will help reproduce the same data.\cr
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LDTestimate}}
}
\examples{
 
#-------Example function usage-----------
n <- 100
L <- 1.5
max.c <- 2.5
pi.r <- 0.5
pi.z <- 0.5
lambda <- 1.33
alpha <- 6.67
beta1 <- 0.29
beta2 <- -0.67

# Generate data from SRD
set.seed(123)
data.A1 <- simLDTdata(n,max.c,pi.r,pi.z,
lambda,alpha,beta1,beta2,L)
data.A2 <- simLDTdata(n,max.c,pi.r,pi.z,
lambda,alpha,beta1,beta2,L)
LDTdata <- cbind(X=c(rep(0,n), rep(1,n)), 
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
