\name{algdiv}
\alias{algdiv}
\title{Compute log(gamma(b)/gamma(a+b)) when b >= 8}
\description{
  Computes \deqn{algdiv(a,b) := \log \frac{\Gamma(b)}{\Gamma(a+b)} = \log
    \Gamma(b) - \log\Gamma(a+b) = \code{lgamma(b) - lgamma(a+b)}}
  in a numerically stable way.

  This is an auxiliary function ins TOMS 708 implementation of
  \code{\link{pbeta}}, aka the incomplete beta function ratio.
}
\usage{
algdiv(a, b)
}
\arguments{
  \item{a, b}{numeric vectors which will be recycled to the same length.}
}
\details{
  Note that this is also useful to compute the Beta function
  \deqn{B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}.}{B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b).}
  Clearly,  \deqn{\log B(a,b) = \log\Gamma(a) + algdiv(a,b) = \log\Gamma(a) - logQab(a,b)}

\preformatted{
 In our ../tests/qbeta-dist.R  we look into computing  log(p*Beta(p,q)) accurately for p << q
        ---------------------
 We are proposing a nice solution there.

 How is this related to algdiv() ?
 }
}
\value{
  a numeric vector of length \code{max(length(a), length(b))} (if neither
  is of length 0, in which case the result has length 0 as well).
}
\references{
  Didonato, A. and Morris, A., Jr, (1992)
  Algorithm 708: Significant digit computation of the incomplete beta
  function ratios,
  \emph{ACM Transactions on Mathematical Software} \bold{18}, 360--373.
}
\author{Didonato, A. and Morris, A., Jr, (1992); \code{algdiv()}'s C
  version from the \R sources, authored by the \R core team; C and \R
  interface: Martin Maechler
}
\seealso{
  \code{\link{gamma}}, \code{\link{beta}};
  my own \code{\link{logQab_asy}()}.
}
\examples{
Qab <- algdiv(2:3, 8:14)
cbind(a = 2:3, b = 8:14, Qab) # recycling with a warning

## algdiv()  and my  logQab_asy()  give *very* similar results for largish b:
all.equal( -   algdiv(3, 100),
           logQab_asy(3, 100), tol=0) # 1.283e-16 !!
(lQab <- logQab_asy(3, 1e10))
## relative error
1 + lQab/ algdiv(3, 1e10) # 0 (64b F 30 Linux; 2019-08-15)
}
\keyword{math}
