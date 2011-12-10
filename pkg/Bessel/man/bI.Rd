\name{bI}
\alias{bI}
\title{Bessel I() function Simple Series Representation}
\description{
  Computes the modified Bessel \eqn{I} function, using one of its basic
  definitions as an infinite series.  The implementation is pure \R,
  working for \code{\link{numeric}}, \code{\link{complex}}, but also
  e.g., for objects of class \code{"\link[Rmpfr:mpfr-class]{mpfr}"}
  from package \pkg{Rmpfr}.
}
\usage{
bI(x, nu, nterm = 800, expon.scaled = FALSE, log = FALSE,
   Ceps = if (isNum) 8e-16 else 2^(-x@.Data[[1]]@prec))
}
\arguments{
  \item{x}{numeric of complex vector, or of another \code{\link{class}}
    for which arithmetic methods are defined, notably
    objects of class \code{\linkS4class{mpfr}}.}
  \item{nu}{non-negative numeric (scalar).}
  \item{nterm}{integer indicating the number of terms to be used.
    should be in the order of \code{abs(x)}, but can be smaller for
    large x.  A warning is given, when \code{nterm} was chosen too
    small.}
  \item{expon.scaled}{logical indicating if the result should be scaled
    by \eqn{exp(-abs(x))}.}
  \item{log}{logical indicating if the logarithm \eqn{log I.()} is
    required.  \emph{is not yet implemented !}%% FIXME
  }
  \item{Ceps}{a relative error tolerance for checking if \code{nterm}
    has been sufficient.  The default is \dQuote{correct} for double
    precision and also for
    multiprecision objects.}
}
% \details{ _____ FIXME
% }
\value{
  a \dQuote{numeric} (or complex or ...) vector of the same class and
  length as \code{x}.
}
\references{
  Abramowitz, M., and Stegun, I. A. (1955, etc).
  \emph{Handbook of mathematical functions}
  (NBS AMS series 55, U.S. Dept. of Commerce).
}
\author{Martin Maechler}
\seealso{
  This package \code{\link{BesselI}}, \pkg{base} \code{\link{besselI}}, etc
}
\examples{
stopifnot(all.equal(bI     (1:10, 1), # R code
                    besselI(1:10, 1)))# internal C code w/ different algorithm
}
\keyword{math}
