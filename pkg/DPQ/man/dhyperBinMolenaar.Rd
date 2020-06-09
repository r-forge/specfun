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
\arguments{
  \item{x}{ .. }
  \item{m}{ .. }
  \item{n}{ .. }
  \item{k}{ .. }
  \item{log}{\code{\link{logical}} indication if the logarithm
    \code{\link{log}(P)} should be returned (instead of \eqn{P}).}
}
\details{
  . . . 
}
\value{
   . . .
}
\references{
 . . .
}
\author{Martin Maechler}
\seealso{
  \R's own \code{\link{dhyper}()} which uses more sophisticated
  computations.
}
\examples{
## The function is simply defined as
function (x, m, n, k, log = FALSE)
  dbinom(x, size = k, prob = hyper2binomP(x, m, n, k), log = log)
}
\keyword{distribution}
