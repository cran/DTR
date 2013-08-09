\name{LDT.mean.estimator}
\alias{LDT.mean.estimator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating LTD restricted mean estimator
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function provides the mean restricted survival estimates and their estimated standard errors for the dynamic treatment regimes (DTRs) based on the data from sequentially randomized designs as proposed in Lunceford, Davidian and Tsiatis (2002). Because different assignments at the first stage are independent to each other,  the function only calculates estimates for one of the assignments at the first stage. For example, if there are two first-stage therapies A1 and A2 available, the function only calculates estimates for one of the two arms (e.g. A1).
}
\usage{
LDT.mean.estimator(data, L=.Machine$double.xmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame (R, Z, U, delta) representing the data for one of the first-stage assignments from sequentially randomized designs. \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored \cr
}
  \item{L}{
%%     ~~Describe \code{t} here~~
restricted survival time. Default is .Machine$double.xmax, which is the largest double value of R. Set L to a numeric number smaller than the maximum follow-up time if restricted follow-up time up to L is considered.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2.
}
\value{
The function returns mean restricted survival estimates and their standard errors for dynamic treatment regimes A1B1 and A1B2.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
\item{MEAN1}{LDT mean restricted survival estimates for A1B1}
\item{MEAN2}{LDT mean restricted survival estimates for A1B2}
\item{MSE1}{estimated standard error for MEAN1}
\item{MSE2}{estimated standard error for MEAN2}
\item{MCOV12}{estimated covariance for MEAN1 and MEAN2}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
}
\note{
%%  ~~further notes~~
The data frame generated from \code{\link{sim.LDT.data}} is the same as the input data frame for \code{\link{LDT.estimator}} and \code{\link{LDT.mean.estimator}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{DTR}}, \code{\link{sim.LDT.data}}, \code{\link{LDT.estimator}}, \code{\link{DTR.Wald.test}}
}
\examples{
require(DTR)
  
#-------Example function usage-----------
n <- 500
L <- 1.5
max.c <- 2.5
pi.r <- 0.5
pi.z <- 0.5
lambda <- 1.33
alpha <- 6.67
beta1 <- 0.29
beta2 <- -0.67

# Generate A1 data from SRD
set.seed(123)  # Set a seed to reproduce the same data every time running the following command
data <- sim.LDT.data(n,max.c,pi.r,pi.z,lambda,alpha,beta1,beta2,L)

# Calculate LTD restricted mean estimates
LDT.mean.estimator(data, L=1.5)
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
