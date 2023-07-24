\name{phyperBinMolenaar}
\title{HyperGeometric Distribution via Molenaar's Binomial Approximation}
\alias{phyperBinMolenaar}% <- to \link{} to
\alias{phyperBinMolenaar.1}
\alias{phyperBinMolenaar.2}
\alias{phyperBinMolenaar.3}
\alias{phyperBinMolenaar.4}
\description{% ../R/hyper-dist.R
  Compute hypergeometric cumulative probabilities via Molenaar's binomial
  approximations.
  The arguments of these functions are \emph{exactly} those of \R's own
  \code{\link{phyper}()}.
}% ~/R/D/r-devel/R/src/library/stats/man/Hypergeometric.Rd
\usage{
phyperBinMolenaar.1(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBinMolenaar.2(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBinMolenaar.3(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBinMolenaar.4(q, m, n, k, lower.tail = TRUE, log.p = FALSE)

phyperBinMolenaar  (q, m, n, k, lower.tail = TRUE, log.p = FALSE) # Deprecated !
}
\arguments{
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in
    \eqn{0,1,\dots, m+n}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\details{
  Molenaar(1970), as cited in Johnson et al (1992), proposed
  \code{phyperBinMolenaar.1()};
  the other three are just using the mathematical symmetries of the
  hyperbolic distribution, swapping \eqn{k} and \eqn{m}, and
  using \code{lower.tail = TRUE} or \code{FALSE}.
}
\value{
 a \code{\link{numeric}} vector, with the length the maximum of the
 lengths of \code{q, m, n, k}.
}
\references{
  % Norman L. Johnson  and  Samuel Kotz  and  Adrienne W. Kemp
  Johnson, N.L., Kotz, S. and Kemp, A.W. (1992)
  Univariate Discrete Distributions, 2nd ed.; Wiley, \doi{10.1002/bimj.4710360207}.\cr
  Chapter 6, mostly Section \emph{5  Approximations and Bounds}, p.256 ff

  Johnson, N.L., Kotz, S. and Kemp, A.W. (2005)
  Univariate Discrete Distributions, 3rd ed.; Wiley; \doi{10.1002/0471715816}.\cr
  Chapter 6, Section \emph{6.5 Approximations and Bounds}, p.268 ff
}
\author{Martin Maechler}
\seealso{
  \code{\link{phyper}}, the hypergeometric distribution, and \R's own
  \dQuote{exact} computation.
  \code{\link{pbinom}}, the binomial distribution functions.

  Our utility \code{\link{phyperAllBin}()}.
}
\examples{
## The first function is simply
function (q, m, n, k, lower.tail = TRUE, log.p = FALSE)
  pbinom(q, size = k, prob = hyper2binomP(q, m, n, k), lower.tail = lower.tail,
        log.p = log.p)
}
\keyword{distribution}
