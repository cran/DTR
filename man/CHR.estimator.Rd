\name{CHR.estimator}
\alias{CHR.estimator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function for calculating cumulative hazard ratio (CHR) estimator
}
\description{
This function computes the estimates for the cumulative hazard ratios (CHRs) between two different dynamic treatment regimes (DTRs) and their variance estimates at given time points as proposed in Tang and Wahed (2013) [Epub ahead of print]. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is implemented in the current version of the function.
}
\usage{
CHR.estimator(fdata, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdata}{
%%     ~~Describe \code{fdata} here~~
a data frame (X, R, Z, U, delta, V) representing the full data from a simple two stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored \cr
V: covariates to be adjusted. The function allows adjustment for one covariate or more than one covariates \cr
}
  \item{t}{
a vector with time points of interest. For example, t=c(1,3,5) for the CHR estimates at 1, 3, and 5 years respectively
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2. Based on four DTRs, six different cumulative hazard ratios (CHRs) are computed: CHR(A1B2 vs. A1B1), CHR(A2B1 vs. A1B1), CHR(A2B2 vs. A1B1), CHR(A2B1 vs. A1B2), CHR(A2B2 vs. A2B1), and CHR(A2B2 vs. A2B1). The natural logarithms of the CHRs are also computed.
}
\value{
The function returns a list for coefficient estimates, CHR estimates and their variance covariance estimates, natural logarithms of CHR estimates and their variance covariance estimates at each given time point. Each list represents one of the time points of interest in order. For example, if we run EST=CHR.estimator(fdata, t=c(1,3,5)), it would return EST[[1]] for estimates at 1 year, EST[[2]] for estimates at 3 years, and EST[[3]] for estimates at 5 years in order. For estimates at each time point, there are a total of six components (COEF, TIME, CHR, VAR[CHR], LOG(CHR), VAR[LOG(CHR)]) for each list:
\item{COEF}{coefficient estimates for the covariates}
\item{TIME}{time of interest}
\item{CHR}{cumulative hazard ratio estimates for comparing A1B2 to A1B1, A2B1 to A1B1, A2B2 to A1B1, A2B1 to A1B2, A2B2 to A1B2, and A2B2 to A2B1}
\item{VAR[CHR]}{estimated variance covariance matrix for CHR}
\item{LOG(CHR)}{natural logarithms of cumulative hazard ratio estimates for comparing A1B2 to A1B1, A2B1 to A1B1, A2B2 to A1B1, A2B1 to A1B2, A2B2 to A1B2, and A2B2 to A2B1}
\item{VAR[LOG(CHR)]}{estimated variance covariance matrix for LOG(CHR)}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{DTR}}, \code{\link{sim.CHR.data}}, \code{\link{updateBeta}}, \code{\link{CHR.Wald.test}}, \code{\link{DTR.CHR.plot}}

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

# Calculate CHR

est <- CHR.estimator(fdata,t=c(1,3,5))
est[[1]] # estimates at 1 year
est[[2]] # estimates at 3 years
est[[3]] # estimates at 5 years
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


