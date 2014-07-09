\name{DTR.object}
\alias{DTR.object}
\title{
Dynamic treatment regime object 
}
\description{
This class of objects is returned by the \code{DTR} class of functions to represent a list of survival estimates and estimated standard errors for dynamic treatment regimes (DTRs). Objects of this class have methods for the functions \code{print}, \code{summary}, and \code{plot}. 
}

\arguments{
\item{DTR}{
dynamic treatment regimes (i.e. A1B1, A1B2, A2B1, and A2B2)
}
\item{records}{
number of observations for each regime
}
\item{events}{
number of events for each regime
}
\item{censorDTR}{
dynamic treatment regimes corresponding to censoring times
}
\item{censortime}{
censoring times
}
\item{censorsurv}{
the survival estimate corresponding to each censoring time (censortime) and dynamic treatment regime (censorDTR)
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
\item{SURV11}{
the survival estimate for A1B1 at each event time
}
\item{SURV12}{
the survival estimate for A1B2 at each event time
}
\item{SURV21}{
the survival estimate for A2B1 at each event time
}
\item{SURV22}{
the survival estimate for A2B2 at each event time
}
\item{SE11}{
the estimated standard error for SURV11 at each event time
}
\item{SE12}{
the estimated standard error for SURV12 at each event time
}
\item{COV1112}{
the estimated covariance between SURV11 and SURV12 at each event time
}
\item{SE21}{
the estimated standard error for SURV21 at each event time
}
\item{SE22}{
the estimated standard error for SURV22 at each event time
}
\item{COV2122}{
the estimated covariance between SURV21 and SURV22 at each event time  
}
}
\seealso{
\code{\link{LDTestimate}},  \code{\link{WRSEestimate}}, \code{\link{print.DTR}}, 
\code{\link{summary.DTR}}, \code{\link{print.summary.DTR}}, \code{\link{plot.DTR}}
}
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
