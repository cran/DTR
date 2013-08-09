\name{DTR-package}
\alias{DTR-package}
\alias{DTR}
\docType{package}
\title{
Estimation and comparison of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials
}
\description{
This is a package for the estimation and comparison of survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials. In a sequentially randomized design, patients are initially randomized to one of the first-stage therapies. Based on their responses to the first-stage therapy, they are then randomized to one of the second-stage therapies. The second-stage therapy could be a rescue therapy if the response is not favorable, or maintenance therapy if favorable response is achieved. There are treatment sequences resulted from such designs: first-stage therapy -> response -> second-stage therapy. The treatment sequences are also referred to as dynamic treatment regimes (DTRs) or adaptive treatment strategies in the literature. \cr
The estimation functions include \code{\link{LDT.estimator}}, \code{\link{LDT.mean.estimator}}, \code{\link{WRSE.estimator}}, and \code{\link{CHR.estimator}}. \cr
The comparisons functions include \code{\link{DTR.Wald.test}}, \code{\link{DTR.Cox.test}}, \code{\link{CHR.Wald.test}}, and \code{\link{DTR.Logrank.test}}. \cr
The plotting functions include \code{\link{DTR.surv.plot}} and \code{\link{DTR.CHR.plot}}. \cr
The functions for data simulation include \code{\link{sim.LDT.data}}, \code{\link{sim.WRSE.data}}, \code{\link{sim.Cox.data}}, \code{\link{sim.CHR.data}}, and \code{\link{sim.Logrank.data}}. 
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2. \cr

The function \code{\link{sim.LDT.data}} generates data sets from sequentially randomized clinical trials as described in the simulation study of Lunceford, Davidian and Tsiatis (2002). \cr
The function \code{\link{LDT.estimator}} computes the estimates of the survival function and their estimated standard errors for DTRs at given time points as proposed in Lunceford, Davidian and Tsiatis (2002) Equation (3) and Equation (10). \cr
The function \code{\link{LDT.mean.estimator}} computes the mean restricted survival estimates and their standard errors for DTRs as proposed in Lunceford, Davidian and Tsiatis (2002). \cr

The function \code{\link{sim.WRSE.data}} generates data sets from sequentially randomized clinical trials as described in the simulation study of Guo and Tsiatis (2005). \cr
The function \code{\link{WRSE.estimator}} computes the weighted risk set estimator (WRSE) of the survival function and their estimated standard errors for DTRs at given time points as proposed in Guo and Tsiatis (2002) Equation (3) and Equation (16). \cr

The function \code{\link{DTR.Wald.test}} compares the survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials based on the LDT estimator proposed in Lunceford, Davidian and Tsiatis (2002) or the WRSE estimator proposed in Guo and Tsiatis (2005) using the Wald-type tests. \cr

The function \code{\link{DTR.surv.plot}} plots the survival functions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials based on the LDT estimator proposed in Lunceford, Davidian and Tsiatis (2002) or the WRSE estimator proposed in Guo and Tsiatis (2005). \cr

The function \code{\link{sim.Cox.data}} generates a data set from sequentially randomized clinical trials as described in the simulation study of Tang and Wahed (2011). \cr
The function \code{\link{DTR.Cox.test}} compares the survival distributions (i.e. hazard functions) of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials after adjustment for covariates as proposed in Tang and Wahed (2011). \cr

The function \code{\link{sim.CHR.data}} generates a data set from sequentially randomized clinical trials as described in the simulation study of Tang and Wahed (2013) [Epub ahead of print]. \cr
The function \code{\link{CHR.estimator}} computes the estimates for the cumulative hazard ratios (CHRs) between two different dynamic treatment regimes (DTRs) and their variance estimates at given time points as proposed in Tang and Wahed (2013) [Epub ahead of print]. \cr
The function \code{\link{CHR.Wald.test}} compares the cumulative hazard functions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials by calculating the natural logarithms of cumulative hazard ratios (CHRs) and performing the Wald-type tests based on natural logarithms of CHRs as proposed in Tang and Wahed (2013) [Epub ahead of print].\cr
The function \code{\link{DTR.CHR.plot}} plots the estimated cumulative hazard ratios (CHR) or the natural logarithms of the estimated CHRs between two different dynamic treatment regimes (DTRs) from sequentially randomized clinical trials using the CHR estimates proposed in Tang and Wahed (2013) [Epub ahead of print]. \cr

The function \code{\link{sim.Logrank.data}} generates a data set from sequentially randomized clinical trials as described in the simulation study of Kidwell and Wahed (2013). \cr
The function \code{\link{DTR.Logrank.test}} compares the survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials using the weighted logrank tests proposed in an unpublished 2005 PhD thesis from North Carolina State University by X. Guo, Feng and Wahed (2008), and Kidwell and Wahed (2013).


\tabular{ll}{
Package: \tab DTR\cr
Type: \tab Package\cr
Version: \tab 1.2 \cr
Date: \tab 2013-08-01\cr
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
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002 \cr
Guo X: Statistical analysis in two-stage randomization designs in clinical trials. PhD thesis, Department of Statistics, North Carolina State University, 2005 \cr
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005 \cr
Feng W, Wahed AS: Supremum weighted log-rank test and sample size for comparing two-stage adaptive treatment strategies. Biometrika 95:695-707, 2008 \cr
Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary variables. Journal of Applied Statistics 38(12):2925-2938, 2011 \cr
Kidwell KM, Wahed AS: Weighted log-rank statistic to compare shared-path adaptive treatment strategies. Biostatistics, 14(2):299-312, 2013 \cr
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print] \cr
}
\seealso{
\code{\link{sim.LDT.data}}, \code{\link{LDT.estimator}}, \code{\link{LDT.mean.estimator}}, \code{\link{sim.WRSE.data}}, \code{\link{WRSE.estimator}}, \code{\link{DTR.Wald.test}}, \code{\link{DTR.surv.plot}}, \code{\link{sim.Cox.data}}, \code{\link{DTR.Cox.test}}, \code{\link{sim.CHR.data}}, \code{\link{CHR.estimator}}, \code{\link{CHR.Wald.test}}, \code{\link{DTR.CHR.plot}}, \code{\link{sim.Logrank.data}}, \code{\link{DTR.Logrank.test}}
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
