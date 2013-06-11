\name{DTR-package}
\alias{DTR-package}
\alias{DTR}
\docType{package}
\title{
Estimation and comparison of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials
}
\description{
This is a package for the estimation and comparison of survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials. In a sequentially randomized design, patients are initially randomized to one of the first-stage therapies. Based on their responses to the first-stage therapy, they are then randomized to one of the second-stage therapies. The second-stage therapy could be a rescue therapy if the response is not favorable, or maintenance therapy if favorable response is achieved. There are treatment sequences resulted from such designs: first-stage therapy -> response -> second-stage therapy. The treatment sequences are also referred to as dynamic treatment regimes (DTRs) or adaptive treatment strategies in the literature. The function \code{\link{sim.LDT.data}} generates data sets from sequentially randomized clinical trials as described in the simulation study of Lunceford, Davidian and Tsiatis (2002). The function \code{\link{LDT.estimator}} computes the estimates of the survival function and their estimated standard errors for DTRs at given time points as proposed in Lunceford, Davidian and Tsiatis (2002) Equation (3) and Equation (10). The function \code{\link{LDT.mean.estimator}} computes the mean restricted survival estimates and their standard errors for DTRs as proposed in Lunceford, Davidian and Tsiatis (2002).
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2.
\tabular{ll}{
Package: \tab DTR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-06-11\cr
License: \tab GPL (>=2)\cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Xinyu Tang\cr
Biostatistics Program, Department of Pediatrics,\cr
University of Arkansas for Medical Sciences\cr
XTang@uams.edu

Maria Melguizo\cr
Biostatistics Program, Department of Pediatrics,\cr
University of Arkansas for Medical Sciences\cr
MSMelguizocastro@uams.edu

}
\references{
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
}
\seealso{
\code{\link{sim.LDT.data}}, \code{\link{LDT.estimator}}, \code{\link{LDT.mean.estimator}}
}
%\examples{
%
%}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
