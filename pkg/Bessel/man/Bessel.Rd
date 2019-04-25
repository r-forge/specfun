\name{Bessel}
\Rdversion{1.1}
\alias{BesselI}
\alias{BesselJ}
\alias{BesselK}
\alias{BesselY}
\title{Bessel Functions of Complex Arguments I(), J(), K(), and Y()}
\description{
  Compute the Bessel functions I(), J(), K(), and Y(), of complex
  arguments \code{z} and real \code{nu},
}
\usage{
BesselI(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0)
BesselJ(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0)
BesselK(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0)
BesselY(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0)
}
\arguments{
  \item{z}{complex or numeric vector.}
  \item{nu}{numeric (scalar).}
  \item{expon.scaled}{logical indicating if the result should be scaled
    by an exponential factor, typically to avoid under- or over-flow.
    See the \sQuote{Details} about the specific scaling.}
  \item{nSeq}{positive integer; if \eqn{> 1}, computes the result for
    a whole \emph{sequence} of \code{nu} values;\cr
    if \code{nu >= 0},\code{nu, nu+1, ..., nu+nSeq-1},\cr
    if \code{nu < 0}, \code{nu, nu-1, ..., nu-nSeq+1}.}
  \item{verbose}{integer defaulting to 0, indicating the level of verbosity
    notably from C code.}
}
\details{
  The case \code{nu < 0} is handled by using simple formula from
  Abramowitz and Stegun,  see details in \code{\link{besselI}()}.

  The scaling activated by \code{expon.scaled = TRUE} depends on the
  function and the scaled versions are  \describe{
    \item{J():}{\code{BesselJ(z, nu, expo=TRUE)}\eqn{:= \exp(-\left|\Im(z)\right|) J_{\nu}(z)}{%
	:= exp(-|Im(z)|) J[nu](z)}}
    \item{Y():}{\code{BesselY(z, nu, expo=TRUE)}\eqn{ := \exp(-\left|\Im(z)\right|) Y_{\nu}(z)}{%
	:= exp(-|Im(z)|) Y[nu](z)}}
    \item{I():}{\code{BesselI(z, nu, expo=TRUE)}\eqn{ := \exp(-\left|\Re(z)\right|) I_{\nu}(z)}{%
	:= exp(-|Re(z)|) I[nu](z)}}
    \item{K():}{\code{BesselK(z, nu, expo=TRUE)}\eqn{ := \exp(z) K_{\nu}(z)}{%
	:= exp(z) K[nu](z)}}
  }
}
\value{
  a complex or numeric vector (or \code{\link{matrix}} with \code{nSeq}
  columns if \code{nSeq > 1})
  of the same length (or \code{\link{nrow}} when \code{nSeq > 1}) and
  \code{\link{mode}} as \code{z}.
}
\references{
  Abramowitz, M., and Stegun, I. A. (1955, etc).
  \emph{Handbook of mathematical functions}
  (NBS AMS series 55, U.S. Dept. of Commerce),
  \url{http://people.math.sfu.ca/~cbm/aands/}

  Wikipedia (20nn). \emph{Bessel Function},
  \url{https://en.wikipedia.org/wiki/Bessel_function}

  %% See also ../vignettes/Bessel.bib
  D. E. Amos (1986)
  Algorithm 644: A portable package for Bessel functions of a complex
  argument and nonnegative order;
  \emph{ACM Trans. Math. Software} \bold{12}, 3, 265--273.

  D. E. Amos (1983)
  \emph{Computation of Bessel Functions of Complex Argument};  Sand83-0083.

  D. E. Amos (1983)
  \emph{Computation of Bessel Functions of Complex Argument and Large Order}; Sand83-0643.

  D. E. Amos (1985)
  \emph{A subroutine package for Bessel functions of a complex
    argument and nonnegative order}; Sand85-1018.

  Olver, F.W.J. (1974).
  \emph{Asymptotics and Special Functions};
  Academic Press, N.Y., p.420
}
\author{
  Donald E. Amos, Sandia National Laboratories, wrote the original
  fortran code.
  Martin Maechler did the translation to C, and partial cleanup
  (replacing \code{goto}'s), in addition to the \R interface.
}
\seealso{
  The base \R functions \code{\link{besselI}()}, \code{\link{besselK}()}, etc.

  The Hankel functions (of first and second kind),
  \eqn{H_{\nu}^{(1)}(z)} and \eqn{H_{\nu}^{(2)}(z)}: \code{\link{Hankel}}.

  The Airy functions \eqn{Ai()} and \eqn{Bi()} and their first
  derivatives, \code{\link{Airy}}.

  For large \code{x} and/or \code{nu} arguments, algorithm AS~644 is not
  good enough, and the results may overflow to \code{Inf} or underflow
  to zero, such that direct computation of \eqn{\log(I_\nu(x))}{log(I[nu](x))} and
  \eqn{\log(K_\nu(x))}{log(K[nu](x))} are desirable.  For this, we provide
  \code{\link{besselI.nuAsym}()}, \code{\link{besselIasym}()} and
  \code{\link{besselK.nuAsym}(*, log= *)}, based on asymptotic expansions.
}
\examples{
%% FIXME: more examples; e.g.  J(i * z, nu) =  c(nu) * I(z, nu)

## For real small arguments, BesselI() gives the same as base::besselI() :
set.seed(47); x <- sort(round(rlnorm(20), 2))
M <- cbind(x, b = besselI(x, 3), B = BesselI(x, 3))
stopifnot(all.equal(M[,"b"], M[,"B"], tol = 2e-15)) # ~4e-16 even
M

## and this is true also for the 'exponentially scaled' version:
Mx <- cbind(x, b = besselI(x, 3, expon.scaled=TRUE),
               B = BesselI(x, 3, expon.scaled=TRUE))
stopifnot(all.equal(Mx[,"b"], Mx[,"B"], tol = 2e-15)) # ~4e-16 even
}
\keyword{math}
