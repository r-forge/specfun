\name{Bern}
\title{Bernoulli Numbers}
\alias{Bern}
\description{%% ../R/hyper-dist.R
  Return the \eqn{n}-th Bernoulli number \eqn{B_n}, (or \eqn{B_n^+}{Bn+},
  see the reference, where \eqn{B_1 = + \frac 1 2}{B_1 = + 1/2}
}
\usage{
Bern(n, verbose = getOption("verbose", FALSE))
}
\arguments{
  \item{n}{integer, \eqn{n \ge 0}{n >= 0}.}
  \item{verbose}{logical indicating if computation should be traced.}
}
%% \details{
%% }
\value{
  The number \eqn{B_n} of type \code{\link{numeric}}.

  A side effect is the \emph{caching} of computed Bernoulli numbers in the
  hidden \code{\link{environment}} \code{.bernoulliEnv}.
}
\references{\url{https://en.wikipedia.org/wiki/Bernoulli_number}
}
\author{Martin Maechler}
\seealso{
  \code{\link[Rmpfr{Bernoulli}} in \CRANpkg{Rmpfr} in arbitrary precision
  via Riemann's \eqn{\zeta}{zeta} function.

  The next version of package \CRANpkg{gmp} should contain
  \code{\link{BernoulliQ}()} which provides exact Bernoulli numbers as big
  rationals (class \code{"bigq"}).
}
\examples{
(B.0.10 <- vapply(0:10, Bern, 1/2))
## [1]  1.00000000 +0.50000000  0.16666667  0.00000000 -0.03333333  0.00000000
## [7]  0.02380952  0.00000000 -0.03333333  0.00000000  0.07575758
}
\keyword{ math }
