\name{DTR-internal}
\title{Internal DTR functions}
\alias{stat_stepribbon}
\alias{StatStepribbon}
\alias{sub.LDTestimate}
\alias{sub.WRSEestimate}
\alias{updateBeta}

\description{
  Internal DTR functions.
}


\usage{
stat_stepribbon(mapping = NULL, data = NULL, geom = "ribbon", position = "identity")
sub.LDTestimate(pdata, t, L)
sub.WRSEestimate(pdata, t)
updateBeta(beta, V, U, delta, w11, w12, w21, w22)
}

\note{
%%  ~~further notes~~
The functions \code{\link{stat_stepribbon}} and \code{\link{StatStepribbon}} were provided by Jeffrey D. Newmiller \email{jdnewmil@gmail.com}, and later revised by Xinyu Tang for newest version of ggplot2. .
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\details{
  These internal \pkg{DTR} functions are not usually called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{DTR} to the next.
  
  
}
\keyword{internal}