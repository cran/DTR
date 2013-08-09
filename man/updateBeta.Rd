\name{updateBeta}
\alias{updateBeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function for updating the coefficient(s) for covariate(s)
}
\description{
The function updates the partial likelihood estimate of the coefficient(s) for covariate(s) (i.e. beta) by solving the pseudo-score equation of the stratified proportional hazards model with dynamic treatment regimes (DTRs) as strata proposed in Tang and Wahed (2013) [Epub ahead of print]. The Newton-Raphson method is used for computation.
}
\usage{
updateBeta(beta, V, U, delta, w11, w12, w21, w22)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
current estimate of the coefficient(s) for covariate(s)
}
  \item{V}{
covariates included in the model. The function allows for one covariate or more than one covariates
}
  \item{U}{
observed survival time, U is death time if delta=1, and U is censoring time if delta=0
}
  \item{delta}{
censoring indicator, delta=1 for died, and delta=0 for censored
}
  \item{w11}{
weights for dynamic treatment regime A1B1
}
  \item{w12}{
weights for dynamic treatment regime A1B2
}
  \item{w21}{
weights for dynamic treatment regime A2B1
}
  \item{w22}{
weights for dynamic treatment regime A2B2
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is for internal use only and is not to be called by the user. \cr

In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2. 
}
\value{
The function returns an updated estimate of the coefficient(s) for covariate(s).
\item{beta_up}{an updated estimate of the coefficient(s) for covariate(s)}
}
\references{
%% ~put references to the literature/web site here ~
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print]
}
\note{
%%  ~~further notes~~
The function \code{\link{CHR.estimator}} calls the function \code{\link{updateBeta}} to obtain the estimates for the coefficient(s). \cr
The function allows one covariate or more than one covariates, but does not allow no adjustment for covariates. \cr
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{DTR}}, \code{\link{sim.CHR.data}}, \code{\link{CHR.estimator}}, \code{\link{CHR.Wald.test}}, \code{\link{DTR.CHR.plot}}

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


