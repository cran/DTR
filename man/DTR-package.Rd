\name{DTR-package}
\alias{DTR-package}
\alias{DTR}
\docType{package}
\title{
Estimation and comparison of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials
}
\description{
This is a package for the estimation and comparison of survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials. In a sequentially randomized design, patients are initially randomized to one of the first-stage therapies. Based on their responses to the first-stage therapy, they are then randomized to one of the second-stage therapies. The second-stage therapy could be a rescue therapy if the response is not favorable, or maintenance therapy if favorable response is achieved. There are treatment sequences resulted from such designs: first-stage therapy -> response -> second-stage therapy. The treatment sequences are also referred to as dynamic treatment regimes (DTRs) or adaptive treatment strategies or treatment policies in the literature. \cr
The estimation functions include \code{\link{LDTestimate}}, \code{\link{WRSEestimate}}, and \code{\link{CHRestimate}}. \cr
The comparisons functions include \code{\link{contrast_wald}}, \code{\link{contrast_chr}}, \code{\link{PHfit}}, \code{\link{contrast_ph}}, and \code{\link{contrast_logrank}}. \cr
The functions for data simulation include \code{\link{simLDTdata}}, \code{\link{simWRSEdata}}, \code{\link{simPHdata}}, \code{\link{simCHRdata}}, and \code{\link{simLRdata}}. 
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four DTRs: A1B1, A1B2, A2B1, and A2B2. \cr

The function \code{\link{simLDTdata}} generates data sets from sequentially randomized clinical trials as described in the simulation study of Lunceford, Davidian and Tsiatis (2002). \cr
The function \code{\link{LDTestimate}} computes the estimates of the survival function and their estimated standard errors for DTRs at observed event times as proposed in Lunceford, Davidian and Tsiatis (2002) Equation (3) and Equation (10). \cr

The function \code{\link{simWRSEdata}} generates data sets from sequentially randomized clinical trials as described in the simulation study of Guo and Tsiatis (2005). \cr
The function \code{\link{WRSEestimate}} computes the weighted risk set estimator (WRSE) of the survival function and their estimated standard errors for DTRs at observed event times as proposed in Guo and Tsiatis (2002) Equation (3) and Equation (16). \cr

The function \code{\link{contrast_wald}} compares the survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials based on the LDT estimator proposed in Lunceford, Davidian and Tsiatis (2002) or the WRSE estimator proposed in Guo and Tsiatis (2005) using the Wald-type tests. \cr

The function \code{\link{simPHdata}} generates a data set from sequentially randomized clinical trials as described in the simulation study of Tang and Wahed (2011). \cr
The function \code{\link{PHfit}} fits a generalized Cox model as proposed in Tang and Wahed (2011). \cr
The function \code{\link{contrast_ph}} compares the survival distributions (i.e. hazard functions) of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials after adjustment for covariates as proposed in Tang and Wahed (2011). \cr

The function \code{\link{simCHRdata}} generates a data set from sequentially randomized clinical trials as described in the simulation study of Tang and Wahed (2013) [Epub ahead of print]. \cr
The function \code{\link{CHRestimate}} computes the estimates for the cumulative hazard ratios (CHRs) between two different dynamic treatment regimes (DTRs) and their variance estimates at given time points as proposed in Tang and Wahed (2013) [Epub ahead of print]. \cr
The function \code{\link{contrast_chr}} compares the cumulative hazard functions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials by calculating the natural logarithms of cumulative hazard ratios (CHRs) and performing the Wald-type tests based on natural logarithms of CHRs as proposed in Tang and Wahed (2013) [Epub ahead of print].\cr

The function \code{\link{simLRdata}} generates a data set from sequentially randomized clinical trials as described in the simulation study of Kidwell and Wahed (2013). \cr
The function \code{\link{contrast_logrank}} compares the survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials using the weighted logrank tests as proposed in Kidwell and Wahed (2013).


\tabular{ll}{
Package: \tab DTR\cr
Type: \tab Package\cr
Version: \tab 1.7 \cr
Date: \tab 2015-12-25\cr
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
Tang X, Melguizo M: DTR: An R Package for Estimation and Comparison of Survival Outcomes of Dynamic Treatment Regimes. Journal of Statistical Software 65(7):1-28 2005 http://www.jstatsoft.org/v65/i07/ \cr
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002 \cr
Guo X: Statistical analysis in two-stage randomization designs in clinical trials. PhD thesis, Department of Statistics, North Carolina State University, 2005 \cr
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005 \cr
Feng W, Wahed AS: Supremum weighted log-rank test and sample size for comparing two-stage adaptive treatment strategies. Biometrika 95:695-707, 2008 \cr
Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary variables. Journal of Applied Statistics 38(12):2925-2938, 2011 \cr
Kidwell KM, Wahed AS: Weighted log-rank statistic to compare shared-path adaptive treatment strategies. Biostatistics, 14(2):299-312, 2013 \cr
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print] \cr
}
\seealso{
\code{\link{simLDTdata}}, \code{\link{LDTestimate}}, \code{\link{simWRSEdata}}, \code{\link{WRSEestimate}}, \code{\link{contrast_wald}}, \code{\link{simPHdata}}, \code{\link{PHfit}}, \code{\link{contrast_ph}}, \code{\link{simCHRdata}}, \code{\link{CHRestimate}}, \code{\link{contrast_chr}}, \code{\link{simLRdata}}, \code{\link{contrast_logrank}}
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
