\name{qntR}
\title{Pure R Implementation of R's qt() / qnt()}
\alias{qntR}
\alias{qntR1}
\description{
  A pure \R implementation of R's C API (\sQuote{Mathlib} specifically)
  \code{qnt()} function which computes (non-central) t quantiles.

  The simple inversion (of \code{pnt()}) scheme has seen to be deficient,
  even in cases where \code{pnt()}, i.e., \R's \code{\link{pt}(.., ncp=*)}
  does not loose accuracy.
}
\usage{
qntR1(p, df, ncp, lower.tail = TRUE, log.p = FALSE,
      pnt = stats::pt, accu = 1e-13, eps = 1e-11)
qntR (p, df, ncp, lower.tail = TRUE, log.p = FALSE,
      pnt = stats::pt, accu = 1e-13, eps = 1e-11)
}
\arguments{
  \item{p, df, ncp}{vectors of probabilities, degrees of freedom, and
    non-centrality parameter; see \code{\link{qt}}.}
  \item{lower.tail, log.p}{logical; see \code{\link{qt}}.}
  \item{pnt}{a \code{\link{function}} for computing the CDF of the (non-central) t-distribution.}
  \item{accu}{a non-negative number, the \dQuote{accu}racy desired in the "root finding" loop.}
  \item{eps}{a non-negative number, used for determining the start interval
    for the root finding.}
}
\value{
  numeric vector of t quantiles, properly recycled in \code{(p, df, ncp)}.
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
  Our \code{\link{qtU}()} and \code{\link{qtAppr}()}; non-central density and probability
  approximations in \code{\link{dntJKBf}}, and e.g., \code{\link{pntR}}.
  Further, \R's \code{\link{qt}}.
}
\examples{
## example where qt() and qntR() "fail" {warnings; --> Inf}
lp <- seq(-30, -24, by=1/4)
summary(p <- exp(lp))
(qp <- qntR( p, df=35, ncp=-7, lower.tail=FALSE))
qp2 <- qntR(lp, df=35, ncp=-7, lower.tail = FALSE, log.p=TRUE)
all.equal(qp, qp2)## same warnings, same values
}
