\name{CHR.Wald.test}
\alias{CHR.Wald.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) based on cumulative hazard ratios (CHRs)
}
\description{
This function compares the cumulative hazard functions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials by calculating the natural logarithms of cumulative hazard ratios (CHRs) and performing the Wald-type tests based on natural logarithms of CHRs as proposed in Tang and Wahed (2013) [Epub ahead of print]. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is implemented in the current version of the function. \cr

}
\usage{
CHR.Wald.test(fdata, t=quantile(fdata$U, 0.75))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdata}{
%%     ~~Describe \code{data} here~~
a data frame (X, R, Z, U, delta, V) representing the full data from a simple two stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored \cr
V: covariates to be adjusted. The function allows adjustment for one covariate or more than one covariates \cr
}
  \item{t}{
%%     ~~Describe \code{t} here~~
a time point of interest. For example, t=5 for the comparisons of survival estimates at 5 years among DTRs. Default t is set to be the 75th percentile of the observed time (i.e. 75th percentile of U).
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2. \cr
Two different comparisons are performed: \cr
1) An overall comparison H0: A1B1=A1B2=A2B1=A2B2 \cr
2) Pairwise comparisons including H0: A1B1=A1B2, H0: A1B1=A2B1, H0: A1B1=A2B2, H0: A1B2=A2B1, H0: A1B2=A2B2, and H0: A2B1=A2B2
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: H0, test statistic, and p-value.
\item{H0}{the null hypotheses being tested, for example, H0 (t=3) for comparisons made at 3 years}
\item{test statistic}{the calculated chi-square test statistic}
\item{p}{the resulting p-value}
}
\references{
%% ~put references to the literature/web site here ~
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print]
}
\note{
%%  ~~further notes~~
The data frame generated from \code{\link{sim.CHR.data}} is the same as the input data frame for \code{\link{CHR.estimator}} and \code{\link{CHR.Wald.test}}. \cr
The function allows one covariate or more than one covariates, but does not allow no adjustment for covariates. \cr
}
\seealso{
\code{\link{DTR}}, \code{\link{sim.CHR.data}}, \code{\link{updateBeta}}, \code{\link{CHR.estimator}}, \code{\link{DTR.CHR.plot}}
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

CHR.Wald.test(fdata, t=3) # Testing
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
