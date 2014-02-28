\name{CHR.object}
\alias{CHR.object}
\title{
Cumulative hazard ratio object 
}
\description{
This class of objects is returned by the \code{CHR} class of functions to represent a list of cumulative hazard ratio estimates and estimated standard errors for dynamic treatment regimes (DTRs). Objects of this class have methods for the functions \code{print}, \code{summary}, and \code{plot}. 
}

\arguments{
\item{coefficients}{
coefficient estimate(s) for the covariate(s)
}
\item{comparison}{
comparisons between DTRs (i.e. A1B2 vs. A1B1, A2B1 vs. A1B1, A2B2 vs. A1B1, A2B1 vs. A1B2, A2B2 vs. A1B2, and A2B2 vs. A2B1)
}
\item{time75P}{
75th percentile of the observed times
}
\item{time}{
event times
}
\item{n.risk}{
number of patients at risk at each event time
}
\item{n.event}{
number of events at each event time
}
\item{CHR1211}{
the cumulative hazard ratio estimate for comparing A1B2 to A1B1 at each event time
}
\item{CHR2111}{
the cumulative hazard ratio estimate for comparing A2B1 to A1B1 at each event time
}
\item{CHR2211}{
the cumulative hazard ratio estimate for comparing A2B2 to A1B1 at each event time
}
\item{CHR2112}{
the cumulative hazard ratio estimate for comparing A2B1 to A1B2 at each event time
}
\item{CHR2212}{
the cumulative hazard ratio estimate for comparing A2B2 to A1B2 at each event time
}
\item{CHR2221}{
the cumulative hazard ratio estimate for comparing A2B2 to A2B1 at each event time
}
\item{SE1211}{
the estimated standard error for CHR1211 at each event time
}
\item{SE2111}{
the estimated standard error for CHR2111 at each event time
}
\item{SE2211}{
the estimated standard error for CHR2211 at each event time
}
\item{SE2112}{
the estimated standard error for CHR2112 at each event time
}
\item{SE2212}{
the estimated standard error for CHR2212 at each event time
}
\item{SE2221}{
the estimated standard error for CHR2221 at each event time
}
\item{COV1211_2111}{
the estimated covariance between CHR1211 and CHR2111 at each event time
}
\item{COV1211_2211}{
the estimated covariance between CHR1211 and CHR2211 at each event time
}
\item{COV1211_2112}{
the estimated covariance between CHR1211 and CHR2112 at each event time
}
\item{COV1211_2212}{
the estimated covariance between CHR1211 and CHR2212 at each event time
}
\item{COV1211_2221}{
the estimated covariance between CHR1211 and CHR2221 at each event time
}
\item{COV2111_2211}{
the estimated covariance between CHR2111 and CHR2211 at each event time
}
\item{COV2111_2112}{
the estimated covariance between CHR2111 and CHR2112 at each event time
}
\item{COV2111_2212}{
the estimated covariance between CHR2111 and CHR2212 at each event time
}
\item{COV2111_2221}{
the estimated covariance between CHR2111 and CHR2221 at each event time
}
\item{COV2211_2112}{
the estimated covariance between CHR2211 and CHR2112 at each event time
}
\item{COV2211_2212}{
the estimated covariance between CHR2211 and CHR2212 at each event time
}
\item{COV2211_2221}{
the estimated covariance between CHR2211 and CHR2221 at each event time
}
\item{COV2112_2212}{
the estimated covariance between CHR2112 and CHR2212 at each event time
}
\item{COV2112_2221}{
the estimated covariance between CHR2112 and CHR2221 at each event time
}
\item{COV2212_2221}{
the estimated covariance between CHR2212 and CHR2221 at each event time
}
\item{CHR1211.LOG}{
the log cumulative hazard ratio estimate for comparing A1B2 to A1B1 at each event time
}
\item{CHR2111.LOG}{
the log cumulative hazard ratio estimate for comparing A2B1 to A1B1 at each event time
}
\item{CHR2211.LOG}{
the log cumulative hazard ratio estimate for comparing A2B2 to A1B1 at each event time
}
\item{CHR2112.LOG}{
the log cumulative hazard ratio estimate for comparing A2B1 to A1B2 at each event time
}
\item{CHR2212.LOG}{
the log cumulative hazard ratio estimate for comparing A2B2 to A1B2 at each event time
}
\item{CHR2221.LOG}{
the log cumulative hazard ratio estimate for comparing A2B2 to A2B1 at each event time
}
\item{SE1211.LOG}{
the estimated standard error for CHR1211.LOG at each event time
}
\item{SE2111.LOG}{
the estimated standard error for CHR2111.LOG at each event time
}
\item{SE2211.LOG}{
the estimated standard error for CHR2211.LOG at each event time
}
\item{SE2112.LOG}{
the estimated standard error for CHR2112.LOG at each event time
}
\item{SE2212.LOG}{
the estimated standard error for CHR2212.LOG at each event time
}
\item{SE2221.LOG}{
the estimated standard error for CHR2221.LOG at each event time
}
\item{COV1211_2111.LOG}{
the estimated covariance between CHR1211.LOG and CHR2111.LOG at each event time
}
\item{COV1211_2211.LOG}{
the estimated covariance between CHR1211.LOG and CHR2211.LOG at each event time
}
\item{COV1211_2112.LOG}{
the estimated covariance between CHR1211.LOG and CHR2112.LOG at each event time
}
\item{COV1211_2212.LOG}{
the estimated covariance between CHR1211.LOG and CHR2212.LOG at each event time
}
\item{COV1211_2221.LOG}{
the estimated covariance between CHR1211.LOG and CHR2221.LOG at each event time
}
\item{COV2111_2211.LOG}{
the estimated covariance between CHR2111.LOG and CHR2211.LOG at each event time
}
\item{COV2111_2112.LOG}{
the estimated covariance between CHR2111.LOG and CHR2112.LOG at each event time
}
\item{COV2111_2212.LOG}{
the estimated covariance between CHR2111.LOG and CHR2212.LOG at each event time
}
\item{COV2111_2221.LOG}{
the estimated covariance between CHR2111.LOG and CHR2221.LOG at each event time
}
\item{COV2211_2112.LOG}{
the estimated covariance between CHR2211.LOG and CHR2112.LOG at each event time
}
\item{COV2211_2212.LOG}{
the estimated covariance between CHR2211.LOG and CHR2212.LOG at each event time
}
\item{COV2211_2221.LOG}{
the estimated covariance between CHR2211.LOG and CHR2221.LOG at each event time
}
\item{COV2112_2212.LOG}{
the estimated covariance between CHR2112.LOG and CHR2212.LOG at each event time
}
\item{COV2112_2221.LOG}{
the estimated covariance between CHR2112.LOG and CHR2221.LOG at each event time
}
\item{COV2212_2221.LOG}{
the estimated covariance between CHR2212.LOG and CHR2221.LOG at each event time
}
}
\seealso{
\code{\link{CHRestimate}},  \code{\link{print.CHR}}, 
\code{\link{summary.CHR}}, \code{\link{print.summary.CHR}}, \code{\link{plot.CHR}}
}
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
