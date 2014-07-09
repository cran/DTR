\name{contrast_logrank}
\alias{contrast_logrank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare dynamic treatment regimes (DTRs) using weighted logrank tests
}
\description{
This function compares the survival distributions of dynamic treatment regimes (DTRs) from sequentially randomized clinical trials using the weighted logrank tests proposed in Kidwell and Wahed (2013). \cr
}
\usage{
contrast_logrank(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame (X, TR, R, Z, U, delta) representing the data from a two-stage randomization design with therapies A1 and A2 available at the first stage, and B1 and B2 available at the second stage. \cr
X: first-stage indicator, X=0 if assigned to A1, and X=1 if assigned to A2 \cr
TR: optional time to response for responders (R=1), only needed if method="WRSE" \cr
R: response status, R=1 for responders, and R=0 for non-responders \cr
Z: second-stage indicator, Z=0 if assigned to B1, and Z=1 if assigned to B2 \cr
U: observed survival time, U is death time if delta=1, and U is censoring time if delta=0 \cr
delta: censoring indicator, delta=1 for died, and delta=0 for censored 
}
}
\details{
In sequentially randomized designs, there could be more than two therapies available at each stage. For simplicity, and to maintain similarity to the most common sequentially randomized clinical trials, a simple two-stage randomization design allowing two treatment options at each stage is used in the current version of the package. In detail, patients are initially randomized to either A1 or A2 at the first stage. Based on their response status, they are then randomized to either B1 or B2 at the second stage. Therefore, there are a total of four dynamic treatment regimes (DTRs): A1B1, A1B2, A2B1, and A2B2. \cr
Two different comparisons are performed: \cr
1) An overall comparison H0: A1B1=A1B2=A2B1=A2B2 \cr
2) Pairwise comparisons including H0: A1B1=A1B2, H0: A1B1=A2B1, H0: A1B1=A2B2, H0: A1B2=A2B1, H0: A1B2=A2B2, and H0: A2B1=A2B2
}
\value{
%%  ~Describe the value returned
The function returns a data frame with four columns: H0, (standardized) test statistic, df, and p.
\item{H0}{the null hypotheses being tested, for example, H0: A1B1=A1B2=A2B1=A2B2 for overall comparison}
\item{(standardized) test statistic}{the calculated (standardized) test statistic, chi-square test statistic for overall comparison, and z test statistic for pairwise comparisons}
\item{df}{the degree of freedom}
\item{p}{the resulting p-value}
}
\references{
%% ~put references to the literature/web site here ~
Guo X: Statistical analysis in two-stage randomization designs in clinical trials. PhD thesis, Department of Statistics, North Carolina State University, 2005 \cr
Feng W, Wahed AS: Supremum weighted log-rank test and sample size for comparing two-stage adaptive treatment strategies. Biometrika 95:695-707, 2008 \cr
Kidwell KM, Wahed AS: Weighted log-rank statistic to compare shared-path adaptive treatment strategies. Biostatistics, 14(2):299-312, 2013 \cr
}
\note{
%%  ~~further notes~~
The data frame generated from \code{\link{simLRdata}} is the same as the input data frame for \code{\link{contrast_logrank}}. \cr
}
\seealso{
\code{\link{simLRdata}}
}
\examples{
\dontrun{
data("LRdata")
contrast_logrank(data=LRdata)}
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
