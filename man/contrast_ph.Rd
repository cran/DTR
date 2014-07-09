\name{contrast_ph}
\alias{contrast_ph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) after adjustment for covariates
}
\description{
This function compares the survival distributions (i.e. hazard functions) of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials after adjustment for covariates as proposed in Tang and Wahed (2011).
}
\usage{
contrast_ph(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fit}{
the result of a call to the \code{\link{PHfit}} function
}
}
\details{
Four different comparisons are performed: \cr
1) An overall comparison H0: A1B1=A1B2=A2B1=A2B2 \cr
2) First-stage comparison H0: A1=A2 \cr
3) Second-stage comparison H0: B1=B2 \cr
4) Pairwise comparisons including H0: A1B1=A1B2, H0: A1B1=A2B1, H0: A1B1=A2B2, H0: A1B2=A2B1, H0: A1B2=A2B2, and H0: A2B1=A2B2
}
\value{
The function returns a data frame with four columns: H0, test statistic, df, and p.
\item{H0}{the null hypotheses being tested, for example, H0 (t=3) if the comparisons are made at t=3}
\item{test statistic}{the calculated chi-square test statistic}
\item{df}{the degree of freedom}
\item{p}{the resulting p-value}
}
\references{
Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary variables. Journal of Applied Statistics 38(12):2925-2938, 2011
}
\seealso{
\code{\link{PHfit}}, \code{\link{coxph.object}}
}
\examples{
\dontrun{
data("PHdata")
f <- PHfit(data=PHdata, covar="V")
contrast_ph(f)}
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
