\name{PHfit}
\alias{PHfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for fitting a generalized proportional hazards model
}
\description{
The function fits a generalized proportional hazards model as proposed in Tang and Wahed (2011). 
}
\usage{
PHfit(data, covar=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame (X, TR, R, Z, U, delta, \ldots) representing the data from a two-stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
TR: time to response \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored \cr
\ldots: other variables \cr
}
  \item{covar}{
covariate(s) to be adjusted in the model. The default (covar=NULL) fits a model without any covariates \cr  
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2. \cr
}
\value{
The function returns an object of class \code{coxph}. See \code{coxph.object} for details. \cr
}
\references{
Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary variables. Journal of Applied Statistics 38(12):2925-2938, 2011
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simPHdata}}, \code{\link{coxph.object}}
}
\examples{
\dontrun{
data("PHdata")
f <- PHfit(data=PHdata, covar="V")
summary(f)
}
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
