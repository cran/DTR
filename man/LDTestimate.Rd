\name{LDTestimate}
\alias{LDTestimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for calculating LTD estimates
}
\description{
The function computes the survival estimates and estimated standard errors for dynamic treatment regimes (DTRs) at the observed event times as proposed in Lunceford, Davidian and Tsiatis (2002) Equation (3) and Equation (10). 
}
\usage{
LDTestimate(data, L = .Machine$double.xmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame (X, R, Z, U, delta) representing the data from a two-stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: treatment arm, X=0 for A1, and X=1 for A2\cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is event time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for event, and delta=0 for censored \cr
%%     ~~Describe \code{data} here~~
}
  \item{L}{
restricted survival time. Default is .Machine$double.xmax, which is the largest double value of R. Set L to a numeric number smaller than the maximum follow-up time if restricted follow-up time up to L is considered.
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2. \cr
}
\value{
The function returns an object of class \code{DTR}. See \code{DTR.object} for details. \cr
}
\references{
%% ~put references to the literature/web site here ~

Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simLDTdata}}, \code{\link{DTR.object}}, \code{\link{print.DTR}}, 
\code{\link{summary.DTR}}, \code{\link{print.summary.DTR}}, \code{\link{plot.DTR}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
data("LDTdata")
est <- LDTestimate(data=LDTdata)
est}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
