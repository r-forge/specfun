\name{ptnc}
\alias{ptnc}
\alias{qtnc}
\title{CDF and Quantiles for Non-Central t distribution}
\description{
  Find the cumulative distribution or quantiles of a non-central t distribution.
}
\usage{
ptnc(t, df=stop("no df arg"), ncp=0, lower.tail = TRUE)
qtnc(p, df=stop("no df arg"), ncp=0, lower.tail = TRUE)
}
\arguments{
  \item{t}{vector of values at which to evaluate the cdf}
  \item{p}{vector of probabilities at which to evaluate the quantiles}
  \item{df}{"degrees of freedom" parameter of the distribution.}
  \item{ncp}{non-centrality parameter of the distribution}
  \item{lower.tail}{logical indicating if lower (or upper) tail
    probabilities are computed or used, respectively.}
}
\value{
  A vector of probabilities (\code{ptnc}) or quantiles (\code{qtnc})
}
\details{
  All parameters can be vectors and will be replicated to the length of
  the longest.

  Fortran code from the DCDFLIB of Barry W. Brown, James Lovato and
  Kathy Russell (originally from \file{ftp://odin.mdacc.tmc.edu/pub/source/dcdflib.f-1.1-tar.gz})
  with original S wrapper by Donald H. MacQueen (\email{macqueen1@llnl.gov})
}
\references{
  Johnson, Kotz,  Balakrishnan (1995)
  \emph{Continuous Univariate Distributions},
  Vol 2, 2nd Edition; Wiley.  Formula on page 532.
}
\examples{
ptnc( seq(-2,2,by=.25) , df=45, ncp=-.5)
qtnc( p=seq(.2,.8,by=.2) , df=11, ncp=1)
}
\keyword{distribution}
