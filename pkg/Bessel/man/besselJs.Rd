\name{besselJs}
\alias{besselJs}
\title{Bessel J() function Simple Series Representation}
\description{
  Computes the modified Bessel \eqn{J} function, using one of its basic
  definitions as an infinite series, e.g. A. & S., p.360, (9.1.10).  The
  implementation is pure \R, working for \code{\link{numeric}},
  \code{\link{complex}}, but also e.g., for objects of class
  \code{"\link[Rmpfr:mpfr-class]{mpfr}"} from package \pkg{Rmpfr}.
}
\usage{
besselJs(x, nu, nterm = 800, log = FALSE,
         Ceps = if (isNum) 8e-16 else 2^(-x@.Data[[1]]@prec))
}
\arguments{
  \item{x}{numeric or complex vector, or of another \code{\link{class}}
    for which arithmetic methods are defined, notably
    objects of class \code{\linkS4class{mpfr}}.}
  \item{nu}{non-negative numeric (scalar).}
  \item{nterm}{integer indicating the number of terms to be used.
    Should be in the order of \code{abs(x)}, but can be smaller for
    large x.  A warning is given, when \code{nterm} was \emph{possibly} too
    small. (Currently, many of these warnings are wrong, as}
  \item{log}{logical indicating if the logarithm \eqn{log J.()} is
    required.}
  \item{Ceps}{a relative error tolerance for checking if \code{nterm}
    has been sufficient.  The default is \dQuote{correct} for double
    precision and also for
    multiprecision objects.}
}
% \details{ _____ FIXME
% }
\value{
  a \dQuote{numeric} (or complex or \code{"\link[Rmpfr:mpfr-class]{mpfr}"})
  vector of the same class and length as \code{x}.
}
\references{
  Abramowitz, M., and Stegun, I. A. (1955, etc).
  \emph{Handbook of mathematical functions}
  (NBS AMS series 55, U.S. Dept. of Commerce).
  \url{http://people.math.sfu.ca/~cbm/aands/page_360.htm}
}
\author{Martin Maechler}
\seealso{
  This package \code{\link{BesselJ}}, \pkg{base} \code{\link{besselJ}}, etc
}
\examples{
stopifnot(all.equal(besselJs(1:10, 1), # our R code
                    besselJ (1:10, 1)))# internal C code w/ different algorithm

## Large 'nu' ...
x <- (0:20)/4
(bx <- besselJ(x, nu=200))# base R's -- gives (mostly wrong) warnings
if(require("Rmpfr")) { ## Use high precision, notably large exponent range, numbers:
  Bx <- besselJs(mpfr(x, 64), nu=200)
  all.equal(Bx, bx, tol = 1e-15)# TRUE -- warnings were mostly wrong; specifically:
  cbind(bx, Bx)
  signif(asNumeric(1 - (bx/Bx)[19:21]), 4) # only [19] had lost accuracy

  ## With*out* mpfr numbers -- using log -- is accurate (here)
  lbx <- besselJs(     x,      nu=200, log=TRUE)
  lBx <- besselJs(mpfr(x, 64), nu=200, log=TRUE)
  cbind(x, lbx, lBx)
  stopifnot(all.equal(asNumeric(log(Bx)), lbx, tol=1e-15),
	    all.equal(lBx, lbx, tol=4e-16))
} # Rmpfr
}
\keyword{math}
