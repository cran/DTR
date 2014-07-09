\name{CHRestimate}
\alias{CHRestimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function for calculating cumulative hazard ratio (CHR) estimates
}
\description{
This function computes the estimates for the cumulative hazard ratios (CHRs) between two different dynamic treatment regimes (DTRs) and their variance estimates at observed event times as proposed in Tang and Wahed (2013) [Epub ahead of print].
}
\usage{
CHRestimate(data, covar=names(data)[!names(data) 
\%in\% c("X", "R", "Z", "U", "delta")])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{fdata} here~~
a data frame (X, R, Z, U, delta, \ldots) representing the data from a two-stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored \cr
\ldots: other variables \cr
}
  \item{covar}{
covariate(s) to be adjusted. The default uses all the variables in the data other than X, R, Z, U and delta \cr  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2. Based on four DTRs, six different cumulative hazard ratios (CHRs) are computed: CHR(A1B2 vs. A1B1), CHR(A2B1 vs. A1B1), CHR(A2B2 vs. A1B1), CHR(A2B1 vs. A1B2), CHR(A2B2 vs. A2B1), and CHR(A2B2 vs. A2B1). The natural logarithms of the CHRs are also computed.
}
\value{
The function returns an object of class \code{CHR}. See \code{CHR.object} for details. \cr
}
\references{
%% ~put references to the literature/web site here ~
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print]
}
\note{
%%  ~~further notes~~
The data frame generated from \code{\link{simCHRdata}} is the same as the input data frame for \code{\link{CHRestimate}}. \cr
The function allows one covariate or more than one covariates, but does not allow no adjustment for covariates. \cr
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{simCHRdata}}, \code{\link{CHR.object}}, \code{\link{print.CHR}}, 
\code{\link{summary.CHR}}, \code{\link{print.summary.CHR}}, \code{\link{plot.CHR}}
}
\examples{
\dontrun{
data("CHRdata")
est <- CHRestimate(data=CHRdata, covar="V1")
est}
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


