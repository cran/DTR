\name{DTR.Wald.test}
\alias{DTR.Wald.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) using Wald-type tests
}
\description{
This function compares the survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials based on the LDT estimator proposed in Lunceford, Davidian and Tsiatis (2002) or the WRSE estimator proposed in Guo and Tsiatis (2005) using the Wald-type tests. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is implemented in the current version of the function. \cr

}
\usage{
DTR.Wald.test(fdata, t=quantile(fdata$U, 0.75), method="LDT", 
rmean=FALSE, L=.Machine$double.xmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdata}{
%%     ~~Describe \code{data} here~~
a data frame (X, TR, R, Z, U, delta) representing the full data from a simple two stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
TR: optional time to response for responders (R=1), only needed if method="WRSE" \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored 
}
  \item{t}{
%%     ~~Describe \code{t} here~~
a time point of interest. For example, t=5 for the comparisons of survival estimates at 5 years among DTRs. Default t is set to be the 75th percentile of the observed time (i.e. 75th percentile of U). The value of t is only used when mean=FALSE.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
the survival estimates to be used for comparisons. The survival estimator proposed in Lunceford, Davidian and Tsiatis (2002) is used when method="LDT" (default), and the weighted risk set estimator (WRSE) proposed in Guo and Tsiatis (2005) is used when method="WRSE".
}
  \item{rmean}{
%%     ~~Describe \code{rmean} here~~
if rmean=FALSE (default), the comparisons are made for the survival estimates at time t. If mean=TRUE, the comparisons are made for the restricted mean survival. mean=TRUE is only allowed when method="LDT".
}
  \item{L}{
%%     ~~Describe \code{L} here~~
restricted survival time. Default is .Machine$double.xmax, which is the largest double value of R. Set L to a numeric number smaller than the maximum follow-up time if restricted follow-up time up to L is considered. The value of L is only used when method="LDT".
}

}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2. \cr
Comparisons can be made among DTRs based on two different estimation methods: \cr
1) "LDT" proposed in Lunceford, Davidian and Tsiatis (2002). Two different comparisons are allowed when \code{method="LDT"}: comparisons of survival estimates at certain time point (e.g. 5 years) if \code{rmean=FALSE} (default), and comparisons of restricted mean survival if \code{rmean=TRUE}. Values less than the maximum follow-up time are allowed for \code{L}. \cr
2) "WRSE" proposed in Guo and Tsiatis (2005). Only comparisons of survival estimates at certain time point (e.g. 5 years) is allowed. Options for \code{rmean} and \code{L} are disabled. \cr
Two different comparisons are performed: \cr
1) An overall comparison H0: A1B1=A1B2=A2B1=A2B2 \cr
2) Pairwise comparisons including H0: A1B1=A1B2, H0: A1B1=A2B1, H0: A1B1=A2B2, H0: A1B2=A2B1, H0: A1B2=A2B2, and H0: A2B1=A2B2
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: H0, test statistic, and p-value.
\item{H0}{the null hypotheses being tested, for example, H0 (t=3) if \code{rmean=FALSE} (default)  and t=3, or H0 (mean) if \code{rmean=TRUE} for the overall comparison}
\item{test statistic}{the calculated chi-square test statistic}
\item{p}{the resulting p-value}
}
\references{
%% ~put references to the literature/web site here ~
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002 \cr
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
}
\note{
%%  ~~further notes~~
The data set generated from \code{\link{sim.WRSE.data}} and \code{\link{sim.LDT.data}} is only the data for one of the first-stage assignments. The input data set for \code{\link{DTR.Wald.test}} is the complete data from a simple two-stage randomization design. See example below for obtaining the complete data using \code{\link{sim.WRSE.data}}. \cr
The input data for \code{method="WRSE"} is different from the input data for \code{method="LDT"}. An additional column \code{TR} for the time to response is needed for \code{method="WRSE"} because of the use of a time-varying weights. \cr
\code{rmean} and \code{L} are only used when \code{method="LDT"}.
}
\seealso{
\code{\link{DTR}}, \code{\link{sim.LDT.data}}, \code{\link{LDT.estimator}}, \code{\link{LDT.mean.estimator}}, \code{\link{sim.WRSE.data}}, \code{\link{WRSE.estimator}}, \code{\link{DTR.surv.plot}}
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

# Testing
DTR.Wald.test(fdata, t=1, method="LDT", rmean=FALSE, L=1.5)
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
