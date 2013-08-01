\name{DTR.Cox.test}
\alias{DTR.Cox.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) after adjustment for covariates
}
\description{
This function compares the survival distributions (i.e. hazard functions) of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials after adjustment for covariates as proposed in Tang and Wahed (2011). For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is implemented in the current version of the function. \cr

}
\usage{
DTR.Cox.test(fdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdata}{
%%     ~~Describe \code{data} here~~
a data frame (X, TR, R, Z, U, delta, V) representing the full data from a simple two stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
TR: time to response \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored \cr
V: covariates to be adjusted. The function allows adjustment for no covariates, one covariate, or more than one covariates \cr
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2. \cr
Four different comparisons are performed: \cr
1) An overall comparison H0: A1B1=A1B2=A2B1=A2B2 \cr
2) First-stage comparison H0: A1=A2 \cr
3) Second-stage comparison H0: B1=B2 \cr
4) Pairwise comparisons including H0: A1B1=A1B2, H0: A1B1=A2B1, H0: A1B1=A2B2, H0: A1B2=A2B1, H0: A1B2=A2B2, and H0: A2B1=A2B2
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: H0, test statistic, and p-value.
\item{H0}{the null hypotheses being tested, for example, H0: A1B1=A1B2=A2B1=A2B2}
\item{test statistic}{the calculated chi-square test statistic}
\item{p}{the resulting p-value}
}
\references{
%% ~put references to the literature/web site here ~
Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary variables. Journal of Applied Statistics 38(12):2925-2938, 2011
}
\note{
%%  ~~further notes~~
The data frame generated from \code{\link{sim.Cox.data}} is the same as the input data frame for \code{\link{DTR.Cox.test}}. \cr

}
\seealso{
\code{\link{DTR}}, \code{\link{sim.Cox.data}}
}
\examples{

require(DTR)
  
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
fdata <- sim.Cox.data(n,max.c,pi.x,pi.z,lambda,alpha,beta1,beta2,beta3,beta4,beta5,gamma)

# Testing
DTR.Cox.test(fdata)
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
