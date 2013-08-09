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

