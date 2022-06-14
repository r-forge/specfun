\name{qtAppr}
\alias{qtAppr}
\title{Compute Approximate Quantiles of Non-Central t Distribution}
\description{
  Compute quantiles (inverse distribution values) for the non-central t distribution.
  using Johnson,Kotz,.. p.521, formula (31.26 a) (31.26 b) & (31.26 c)

  Note that \code{\link{qt}(.., ncp=*)} did not exist yet in 1999, when MM
  implemented \code{qtAppr()}.
}
\usage{
qtAppr(p, df, ncp, lower.tail = TRUE, log.p = FALSE, method = c("a", "b", "c"))
}
\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/TDist.Rd : qt(p, df, ncp, lower.tail = TRUE, log.p = FALSE)
  \item{p}{vector of probabilities.}
  \item{df}{degrees of freedom \eqn{> 0}, maybe non-integer.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qt}()}.}
  \item{method}{a string specifying the approximation method to be used.}
}
%% \details{
%% }
\value{
  ...
}
\references{
  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995)
  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley.
  \cr Chapter 31, Section \emph{6  Approximation}, p.519 ff
}
\author{Martin Maechler,  6 Feb 1999}
\seealso{
  Our \code{\link{qtU}()}; several non-central density and probability
  approximations in \code{\link{dntJKBf}}, and e.g., \code{\link{pntR}}.
  Further, \R's \code{\link{qt}}.
}
\examples{
  ## TODO -- currently see ../tests/t-nonc-tst.R
}
\keyword{distribution}
