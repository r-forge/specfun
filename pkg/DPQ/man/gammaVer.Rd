\name{gammaVer}
\alias{gammaVer}
\title{Gamma Function Versions}
\description{
  Provide different variants or versions of computing the Gamma
  (\eqn{\Gamma}) function.
}
\usage{% >>> ../R/qgamma-fn.R  ----  ../src/gamma-variants.c <<<------
gammaVer(x, version, traceLev = 0)
}
\arguments{
  \item{x}{numeric vector of absissa value for the Gamma function.}
  %% FIXME: use c("..", ...) match.arg() etc:
  \item{version}{integer in {1,2,..,5} specifying which variant is desired.}
  \item{traceLev}{non-negative integer indicating the amount of diagnostic
    \dQuote{tracing} output to the console during computation.}
}
\details{
  All of these are good algorithms to compute \eqn{\Gamma(x)} (for real
  \eqn{x}), and indeed correspond to the versions \R's implementation of
  \code{\link{gamma}(x)} over time.  More specifically, the current version
  numbers correspond to
  \enumerate{
    \item .
    \item .
    \item .
    \item Used in \R from ... up to versions 4.2.z
    \item Possibly to be used in \R 4.3.z and newer.
  }
}
\value{
  numeric vector as \code{x}
}
\references{
 .... TODO ....
}
\author{Martin Maechler}
\seealso{
 \code{\link{gamma}()}, \R's own Gamma function.
}
\examples{
xx <- seq(-4, 10, by=1/2)
gx <- sapply(1:5, gammaVer, x=xx)
gamx <- gamma(xx)
cbind(xx, gx, gamma=gamx)
apply(gx, 2, all.equal, target=gamx, tol = 0) # typically: {T,T,T,T, 1.357e-16}
stopifnot( apply(gx, 2, all.equal, target = gamx, tol = 1e-14))
                                                 # even 2e-16 (Lnx, 64b, R 4.2.1)
}
\keyword{math}
