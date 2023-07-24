\name{dhyperBinMolenaar}
\alias{dhyperBinMolenaar}
\title{HyperGeometric (Point) Probabilities via Molenaar's Binomial Approximation}
\description{
  Compute hypergeometric (point) probabilities via Molenaar's binomial
  approximation, \code{\link{hyper2binomP}()}.
}
\usage{
dhyperBinMolenaar(x, m, n, k, log = FALSE)
}
\arguments{% see  ~/R/D/r-devel/R/src/library/stats/man/Hypergeometric.Rd
  % exactly the same arguments as those in \code{\link{dhyper}().}
  \item{x}{(vector of) the number of white balls drawn without replacement
    from an urn which contains both black and white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence in \eqn{0,1,\dots, m+n}.}
  \item{log}{\code{\link{logical}} indication if the logarithm
    \code{\link{log}(P)} should be returned (instead of \eqn{P}).}
}
\value{
  a \code{\link{numeric}} vector, with the length the maximum of the
  lengths of \code{x, m, n, k}.
}
\references{See those in \code{\link{phyperBinMolenaar}}.}
\author{Martin Maechler}
\seealso{
  \code{\link{hyper2binomP}()};
  \R's own \code{\link{dhyper}()} which uses more sophisticated
  computations.
}
\examples{
## The function is simply defined as
function (x, m, n, k, log = FALSE)
  dbinom(x, size = k, prob = hyper2binomP(x, m, n, k), log = log)
}
\keyword{distribution}
