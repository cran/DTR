\name{DTR.Logrank.test}
\alias{DTR.Logrank.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) using weighted logrank tests
}
\description{
This function compares the survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials using the weighted logrank tests proposed in an unpublished 2005 PhD thesis from North Carolina State University by X. Guo, Feng and Wahed (2008), and Kidwell and Wahed (2013). For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is implemented in the current version of the function. \cr

}
\usage{
DTR.Logrank.test(fdata)
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
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2. \cr
Two different comparisons are performed: \cr
1) An overall comparison H0: A1B1=A1B2=A2B1=A2B2 \cr
2) Pairwise comparisons including H0: A1B1=A1B2, H0: A1B1=A2B1, H0: A1B1=A2B2, H0: A1B2=A2B1, H0: A1B2=A2B2, and H0: A2B1=A2B2
}
\value{
%%  ~Describe the value returned
The function returns a data set with columns: H0, (standardized) test statistic, and p-value.
\item{H0}{the null hypotheses being tested, for example, H0: A1B1=A1B2=A2B1=A2B2 for overall comparison}
\item{(standardized) test statistic}{the calculated (standardized) test statistic, chi-square test statistic for overall comparison, and z test statistic for pairwise comparisons}
\item{p}{the resulting p-value}
}
\references{
%% ~put references to the literature/web site here ~
Guo X: Statistical analysis in two-stage randomization designs in clinical trials. PhD thesis, Department of Statistics, North Carolina State University, 2005 \cr
Feng W, Wahed AS: Supremum weighted log-rank test and sample size for comparing two-stage adaptive treatment strategies. Biometrika 95:695-707, 2008 \cr
Kidwell KM, Wahed AS: Weighted log-rank statistic to compare shared-path adaptive treatment strategies. Biostatistics, 14(2):299-312, 2013 \cr
}
\note{
%%  ~~further notes~~
The data frame generated from \code{\link{sim.Logrank.data}} is the same as the input data frame for \code{\link{DTR.Logrank.test}}. \cr
}
\seealso{
\code{\link{DTR}}, \code{\link{sim.Logrank.data}}
}
\examples{

require(DTR)
  
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
fdata <- sim.Logrank.data(n,max.c,pi.x,pi.r,pi.z,mean.NR.1,mean.NR.2,
mean.R.1,mean.R.2,mean.RE.11,mean.RE.12,mean.RE.21,mean.RE.22)

# Testing
DTR.Logrank.test(fdata)
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
