\name{besselI.nuAsym}
\title{Asymptotic Expansion of Bessel I(x,nu) and K(x,nu) for Large nu (and x)}
\alias{besselI.nuAsym}
\alias{besselK.nuAsym}
\description{
  Compute Bessel functions \eqn{I_{\nu}(x)}{I[nu](x)} and
  \eqn{K_{\nu}(x)}{K[nu](x)} for large \eqn{\nu}{nu} and possibly large
  \eqn{x},
  using asymptotic expansions in Debye polynomials.
}
\usage{
besselI.nuAsym(x, nu, k.max, expon.scaled = FALSE, log = FALSE)
besselK.nuAsym(x, nu, k.max, expon.scaled = FALSE, log = FALSE)
}
\arguments{
  \item{x}{numeric or \code{\link{complex}}, with real part \eqn{\ge 0}{>= 0}.}
  \item{nu}{numeric; The \emph{order} (maybe fractional!) of the
    corresponding Bessel function.}
  \item{k.max}{integer number of terms in the expansion.  Must be in
    \code{0:5}, currently.}
  \item{expon.scaled}{logical; if \code{TRUE}, the results are
    exponentially scaled, the same as in the corresponding
    \code{\link{BesselI}()} and \code{BesselK()} functions in order to avoid overflow
    (\eqn{I_{\nu}}{I(nu)}) or underflow (\eqn{K_{\nu}}{K(nu)}), respectively.}
  \item{log}{logical; if TRUE, \eqn{\log(f(.))}{log(f(.))} is returned instead of \eqn{f}.}
}
\details{
  Abramowitz & Stegun , page 378, has formula 9.7.7 and 9.7.8 for the asymptotic
  expansions of \eqn{I_{\nu}(x)} and \eqn{K_{\nu}(x)}, respectively,
  also saying \emph{
    When \eqn{\nu \to +\infty}{nu -> oo}, these expansions
    (of \eqn{I_{\nu}(\nu z)}{I[nu](nu * z)} and
    \eqn{K_{\nu}(\nu z)}{K[nu](nu * z)})
    hold uniformly with respect to \eqn{z} in the sector
    \eqn{|arg z| \le \frac{1}{2} \pi - \epsilon}{|arg z| <= pi/2 - eps},
    where \eqn{\epsilon}{eps} iw qn arbitrary positive number.}
  and for this reason, we require \eqn{\Re(x) \ge 0}{Re(x) >= 0}.

  The Debye polynomials \eqn{u_k(x)} are defined in 9.3.9 and 9.3.10 (page 366).
}
\value{
  a numeric vector of the same length as the long of \code{x} and
  \code{nu}. (usual argument recycling is applied implicitly.)
}
\references{
  Abramowitz, M., and Stegun, I. A. (1964, etc).
  \emph{Handbook of mathematical functions}, pp. 366, 378.
}
\author{Martin Maechler}
\seealso{
  From this package \pkg{Bessel}: \code{\link{BesselI}()}; further,
  \code{\link{besselIasym}()} for the case when \eqn{x} is large and
  \eqn{\nu} is small or moderate.

  Further, from \pkg{base}: \code{\link{besselI}}, etc.
}
\examples{
x <- c(1:10, 20, 50, 100, 100000)
nu <- c(1, 10, 20, 50, 10^(2:10))

%% FIXME: do something more sensible -- compare with  besselI() and *K()
sapply(0:4, function(k.)
            sapply(nu, function(n.)
                   besselI.nuAsym(x, nu=n., k.max = k., log = TRUE)))

sapply(0:4, function(k.)
            sapply(nu, function(n.)
                   besselK.nuAsym(x, nu=n., k.max = k., log = TRUE)))
}
\keyword{math}
