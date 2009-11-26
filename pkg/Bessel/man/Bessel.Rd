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
BesselI(z, nu, expon.scaled = FALSE, nSeq = 1)
BesselJ(z, nu, expon.scaled = FALSE, nSeq = 1)
BesselK(z, nu, expon.scaled = FALSE, nSeq = 1)
BesselY(z, nu, expon.scaled = FALSE, nSeq = 1)
}
\arguments{
  \item{z}{complex or numeric vector.}
  \item{nu}{numeric (scalar).}
  \item{expon.scaled}{logical indicating if the result should be scaled
    by an exponential factor (typically to avoid under- or over-flow).}
  \item{nSeq}{positive integer; if \eqn{> 1}, computes the result for
    a whole \emph{sequence} of \code{nu} values;\cr
    if \code{nu >= 0},\code{nu, nu+1, ..., nu+nSeq-1},\cr
    if \code{nu < 0}, \code{nu, nu-1, ..., nu-nSeq+1}.}
}
\details{
%% FIXME: explain  expon.scaled ---> ../src/zbsubs.f

The case \code{nu < 0} is handled by using simple formula from
Abramowitz and Stegun. %% FIXME : be explicit
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
  (NBS AMS series 55, U.S. Dept. of Commerce).

  D. E. Amos (1986)
  A portable package for Bessel functions of a complex argument
  and nonnegative order;
  \emph{ACM Trans. Math. Software} \bold{12}, 3, 265--273.

  D. E. Amos (1983)
  \emph{Computation of Bessel Functions of Complex Argument};
  Sand83-0083.

  D. E. Amos (1983)
  \emph{Computation of Bessel Functions of Complex Argument and
    Large Order}; Sand83-0643.

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
  Martin Maechler did the \R interface.
}
\seealso{
  The base \R functions \code{\link{besselI}}, etc.
}
\examples{
%% FIXME: more examples; e.g.  J(i * z, nu) =  c(nu) * I(z, nu)

## For real small arguments, BesselI() gives the same as base::besselI() :
set.seed(47); x <- sort(round(rlnorm(20), 2))
M <- cbind(x, b = besselI(x, 3), B = BesselI(x, 3))
stopifnot(all.equal(M[,"b"], M[,"B"]))
M
}
\keyword{math}
