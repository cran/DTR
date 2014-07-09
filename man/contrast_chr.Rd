\name{contrast_chr}
\alias{contrast_chr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) based on cumulative hazard ratios (CHRs)
}
\description{
This function compares the cumulative hazard functions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials by calculating the natural logarithms of cumulative hazard ratios (CHRs) and performing the Wald-type tests based on natural logarithms of CHRs as proposed in Tang and Wahed (2013) [Epub ahead of print].
}
\usage{
contrast_chr(est, t = quantile(est$time, 0.75))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{est}{
the result of a call to the \code{\link{CHRestimate}} function
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
Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print]
}
\seealso{
\code{\link{CHRestimate}},  \code{\link{CHR.object}}
}
\examples{
\dontrun{
data("CHRdata")
est <- CHRestimate(data=CHRdata, covar="V1")
contrast_chr(est, t=1.5)}
}
\keyword{survival analysis}
\keyword{sequentially randomized clinical trial}
\keyword{sequentially randomized design}
\keyword{treatment sequence}
\keyword{dynamic treatment regime}
\keyword{adaptive treatment strategy}
\keyword{Inverse weighting}
