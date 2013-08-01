\name{WRSE.estimator}
\alias{WRSE.estimator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function for calculating WRSE estimator
}
\description{
The function computes the weighted risk set estimates (WRSE) of the survival functions and their estimated standard errors for dynamic treatment regimes (DTRs) at given time points as proposed in Guo and Tsiatis (2005) Equation (3) and Equation (16). Because different assignments at the first stage are independent to each other,  the function only calculates estimates for one of the assignments at the first stage. For example, if there are two first-stage therapies A1 and A2 available, the function only calculates estimates for one of the two arms (e.g. A1). 


}
\usage{
WRSE.estimator(data, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame (TR, R, Z, U, delta) representing the data for one of the first-stage assignments from sequentially randomized designs. \cr
TR: the time to response for responders (R=1) \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0  \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored 

}
  \item{t}{
a vector with time points of interest. For example, t=c(150,500,700) for the survival estimates at 150, 500, and 700 days respectively
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2.
}
\value{
The function returns a matrix for the WRSEs and their standard errors for dynamic treatment regimes A1B1 and A1B2 at given time points. Each row represents one of the time points of interest in order. There are a total of six columns (TIME, SURV1, SURV2, SE1, SE2, COV12).
\item{TIME}{times of interest}
\item{SURV1}{WRSE for A1B1}
\item{SURV2}{WRSE for A1B2}
\item{SE1}{estimated standard error for SURV1}
\item{SE2}{estimated standard error for SURV2}
\item{COV12}{estimated covariance for SURV1 and SURV2}

}
\references{
%% ~put references to the literature/web site here ~
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
}
\note{
%%  ~~further notes~~
The data frame generated from \code{\link{sim.WRSE.data}} is the same as the input data frame for \code{\link{WRSE.estimator}}. \cr
If the values for \code{max.c}, \code{mean.T0}, \code{mean.TR}, \code{mean.T1} and \code{mean.T2} were given in days, in months, or in years when simulating the data, the values for the \code{t} should be given in days, in months, or in years accordingly.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{DTR}}, \code{\link{sim.WRSE.data}}, \code{\link{DTR.Wald.test}}, \code{\link{DTR.surv.plot}}

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

# Generate data as in Guo and Tsiatis (2005)
set.seed(123)  # Set a seed to reproduce the same data every time running the following command
data <- sim.WRSE.data(n,max.c,pi.r,pi.z,mean.T0,mean.TR,mean.T1,mean.T2)

# Calculate WRSE
WRSE.estimator(data,t=c(150,500,700))
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


