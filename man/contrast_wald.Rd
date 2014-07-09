\name{contrast_wald}
\alias{contrast_wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) using Wald-type tests
}
\description{
This function compares the survival estimates at specific time point among dynamic treatment regimes (DTRs) using the Wald-type tests.
}
\usage{
contrast_wald(est, t = quantile(est$time, 0.75))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{est}{
the result of a call to the \code{\link{LDTestimate}} function or \code{\link{WRSEestimate}} function
}
\item{t}{
a time point of interest. For example, t=5 for the comparisons of survival estimates at 5 years among DTRs. Default t is set to be the 75th percentile of the observed time (i.e. 75th percentile of U).
}
}
\details{
Two different comparisons are performed: \cr
1) An overall comparison H0: A1B1=A1B2=A2B1=A2B2 \cr
2) Pairwise comparisons including H0: A1B1=A1B2, H0: A1B1=A2B1, H0: A1B1=A2B2, H0: A1B2=A2B1, H0: A1B2=A2B2, and H0: A2B1=A2B2
}
\value{
The function returns a data frame with four columns: H0, test statistic, df, and p.
\item{H0}{the null hypotheses being tested, for example, H0 (t=3) if the comparisons are made at t=3}
\item{test statistic}{the calculated chi-square test statistic}
\item{df}{the degree of freedom}
\item{p}{the resulting p-value}
}
\references{
Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002 \cr
Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
}
\seealso{
\code{\link{LDTestimate}},  \code{\link{WRSEestimate}}, \code{\link{DTR.object}}
}
\examples{
\dontrun{
# LDT estimates
data("LDTdata")
est <- LDTestimate(data=LDTdata)
contrast_wald(est, t=1)

# WRSE estimates
data("WRSEdata")
est <- WRSEestimate(data=WRSEdata)
contrast_wald(est, t=300)
}
}
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
