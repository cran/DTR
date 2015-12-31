\name{NEWS}
\title{NEWS file for the DTR package}
\section{Changes in version 1.2}{
  \itemize{
    \item sim.LDT.data.R: change L to be optional input with default as .Machine$double.xmax, which is the largest double value of R
    \item LDT.estimator.R: change L to be optional input with default as .Machine$double.xmax, which is the largest double value of R
    \item LDT.estimator.R: add a covariance estimate COV12 between SURV1 and SURV2 in the returned value
    \item LDT.estimator.R: change the returned value as a data frame
    \item LDT.mean.estimator.R: change L to be optional input with default as .Machine$double.xmax, which is the largest double value of R
    \item LDT.mean.estimator.R: add a covariance estimate MCOV12 between MEAN1 and MEAN2 in the returned value
    \item LDT.mean.estimator.R: change the returned value as a data frame
    \item New function sim.WRSE.data.R added to simulate data from sequentially randomized clinical trials (Guo and Tsiatis 2005)
    \item New function WRSE.estimator.R added to calculate WRSE estimator (Guo and Tsiatis 2005)
    \item New function DTR.Wald.test.R added to compare dynamic treatment regimes (DTRs) using Wald-type tests
    \item New function DTR.surv.plot.R added to plot survival curves for dynamic treatment regimes (DTRs)
    \item New function sim.Cox.data.R added to simulate data from sequentially randomized clinical trials (Tang and Wahed 2011)
    \item New function DTR.Cox.test.R added to compare dynamic treatment regimes (DTRs) after adjustment for covariates (Tang and Wahed 2011)
    \item New function sim.CHR.data.R added to simulate data from sequentially randomized clinical trials (Tang and Wahed 2013 [Epub ahead of print])
    \item New function updateBeta.R added to update the coefficient estimate(s) for covariate(s)
    \item New function CHR.estimator.R added to calculating cumulative hazard ratio (CHR) estimator (Tang and Wahed 2013 [Epub ahead of print])
    \item New function CHR.Wald.test.R added to compare dynamic treatment regimes (DTRs) based on cumulative hazard ratios (CHRs) (Tang and Wahed 2013 [Epub ahead of print])
    \item New function DTR.CHR.plot.R added to plot (log) cumulative hazard ratio between different dynamic treatment regimes (DTRs)
    \item New function sim.Logrank.data.R added to simulate data from sequentially randomized clinical trials (Kidwell and Wahed 2013)
    \item New function DTR.Logrank.test.R added to compare dynamic treatment regimes (DTRs) using weighted logrank tests
}
}
\section{Changes in version 1.3}{
  \itemize{
    \item The function sim.LDT.data was renamed to simLDTdata  
    \item The function LDT.estimator.R was replaced by the function LDTestimate.R
    \item The function LDTestimate.R returns an object of class DTR containing all the survival estimates and estimated standard errors at observed event times
    \item The function LDT.mean.estimator.R was deleted
    \item The function sim.WRSE.data was renamed to simWRSEdata
    \item The function WRSE.estimator.R was replaced by the function WRSEestimate.R
    \item The function WRSEestimate.R returns an object of class DTR containing all the survival estimates and estimated standard errors at observed event times    
    \item The S3 functions print.DTR, summary.DTR, print.summary.DTR, and plot.DTR were added
    \item The function DTR.Wald.test.R was replaced by the function contrast.wald.R
    \item The function contrast.wald.R implements the function wald.test() from the package aod to perform wald-type tests
    \item The function sim.CHR.data was renamed to simCHRdata
    \item The function CHR.estimator.R was replaced by the function CHRestimate.R
    \item The function CHRestimate.R returns an object of class CHR containing all the cumulative hazard ratio estimates and estimated standard errors at observed event times
    \item The S3 functions print.CHR, summary.CHR, print.summary.CHR, and plot.CHR were added
    \item The function CHR.Wald.test.R was replaced by the function contrast.chr.R 
    \item The function contrast.chr.R implements the function wald.test() from the package aod to perform wald-type tests using the CHR estimates and their variance/covariance estimates
    \item The function sim.Cox.data was renamed to simPHdata
    \item The function DTR.Cox.test was replaced by two separate functions PHfit.R and contrast.ph.R
    \item The function PHfit.R returns an object of class coxph containing the proportional hazard model fit
    \item The function contrast.ph.R implements the function wald.test() from the package aod to perform wald-type tests based on the returned coefficient and variance/covariance estimated from PHfit()
    \item The function sim.Logrank.data.R was renamed to simLRdata
    \item The function DTR.Logrank.test.R was replaced by the function contrast.logrank.R
    \item The data sets LDTdata, WRSEdata, CHRdata, PHdata, and LRdata were added to data folder
    \item datalist was added to data folder
}
}
\section{Changes in version 1.4}{
  \itemize{
    \item The function contrast.wald was renamed to contrast_wald  
    \item The function contrast.chr was renamed to contrast_chr
    \item The function contrast.logrank was renamed to contrast_logrank
    \item The function contrast.ph was renamed to contrast_ph
    \item An input option covar was added to the function PHfit to allow user to specify the covariates to be adjusted
    \item An input option covar was added to the function CHRestimate to allow user to specify the covariates to be adjusted
    \item An input option censored was added to the function plot.DTR for plotting censoring ticks
}
}
\section{Changes in version 1.5}{
  \itemize{
    \item A bug in the function CHRestimate was fixed  
}
}
\section{Changes in version 1.6}{
  \itemize{
    \item Add CITATION file from JSS  
    \item Add reference to DTR-package.RD
}
}
\section{Changes in version 1.7}{
  \itemize{
    \item Delete internal function stairstepn
    \item Internal functions stat_stepribbon and StatStepribbon were revised for newest version of ggplot2
}
}