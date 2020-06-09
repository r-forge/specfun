\name{phyperBinMolenaar}
\title{HyperGeometric Distribution via Molenaar's Binomial Approximation}
\alias{phyperBinMolenaar}
\alias{phyperBinMolenaar.1}
\alias{phyperBinMolenaar.2}
\alias{phyperBinMolenaar.3}
\alias{phyperBinMolenaar.4}
\description{% ../R/hyper-dist.R
  Compute hypergeometric cumulative probabilities via Molenaar's binomial
  approximations.
  The arguments of these functions are \emph{exactly} those of \R's own
  \code{\link{phyper}()}.

  . . .
}% ~/R/D/r-devel/R/src/library/stats/man/Hypergeometric.Rd
\usage{
phyperBinMolenaar  (q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBinMolenaar.1(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBinMolenaar.2(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBinMolenaar.3(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBinMolenaar.4(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
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
%% \details{
%%   . . .
%% }
\value{
 . . .  %% FIXME
}
\references{
  % Norman L. Johnson  and  Samuel Kotz  and  Adrienne W. Kemp
  Johnson, N.L., Kotz, S. and Kemp, A.W. (1992)
  Univariate Discrete Distributions, 2nd ed.; Wiley.
  \cr Chapter 6, mostly Section \emph{5  Approximations and Bounds}, p.256 ff
}
\author{Martin Maechler}
\seealso{
  \code{\link{phyper}}, the hypergeometric distribution, and \R's own
  \dQuote{exact} computation.
  \code{\link{pbinom}}, the binomial distribution functions.
}
\examples{
## The function is currently defined as
function (q, m, n, k, lower.tail = TRUE, log.p = FALSE)
pbinom(q, size = k, prob = hyper2binomP(q, m, n, k), lower.tail = lower.tail,
    log.p = log.p)
}
\keyword{distribution}
