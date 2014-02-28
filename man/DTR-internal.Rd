\name{DTR-internal}
\title{Internal DTR functions}
\alias{stairstepn}
\alias{stat_stepribbon}
\alias{StatStepribbon}
\alias{sub.LDTestimate}
\alias{sub.WRSEestimate}
\alias{updateBeta}

\description{
  Internal DTR functions.
}


\usage{
stairstepn(data, direction = "hv", yvars = "y")
stat_stepribbon(mapping = NULL, data = NULL, geom = "ribbon", position = "identity")
sub.LDTestimate(pdata, t, L)
sub.WRSEestimate(pdata, t)
updateBeta(beta, V, U, delta, w11, w12, w21, w22)
}

\note{
%%  ~~further notes~~
The functions \code{\link{stairstepn}}, \code{\link{stat_stepribbon}}, and \code{\link{StatStepribbon}} were provided by Jeffrey D. Newmiller 
 \email{jdnewmil@gmail.com}.
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\details{
  These internal \pkg{DTR} functions are not usually called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{DTR} to the next.
  
  
}
\keyword{internal}