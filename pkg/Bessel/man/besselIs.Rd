\name{bI}
\alias{bI}
\alias{besselIs}
\title{Bessel I() function Simple Series Representation}
\description{
  Computes the modified Bessel \eqn{I} function, using one of its basic
  definitions as an infinite series.  The implementation is pure \R,
  working for \code{\link{numeric}}, \code{\link{complex}}, but also
  e.g., for objects of class \code{"\link[Rmpfr:mpfr-class]{mpfr}"}
  from package \pkg{Rmpfr}.
}
\usage{
besselIs(x, nu, nterm = 800, expon.scaled = FALSE, log = FALSE,
         Ceps = if (isNum) 8e-16 else 2^(-x@.Data[[1]]@prec))
}
\arguments{
  \item{x}{numeric or complex vector, or of another \code{\link{class}}
    for which arithmetic methods are defined, notably
    objects of class \code{\linkS4class{mpfr}} (package \CRANpkg{Rmpfr}).}
  \item{nu}{non-negative numeric (scalar).}
  \item{nterm}{integer indicating the number of terms to be used.
    Should be in the order of \code{abs(x)}, but can be smaller for
    large x.  A warning is given, when \code{nterm} was chosen too
    small.}
  \item{expon.scaled}{logical indicating if the result should be scaled
    by \eqn{exp(-abs(x))}.}
  \item{log}{logical indicating if the logarithm \eqn{log I.()} is
    required.  This allows even more precision than
    \code{expon.scaled=TRUE} in some cases.}
  \item{Ceps}{a relative error tolerance for checking if \code{nterm}
    has been sufficient.  The default is \dQuote{correct} for double
    precision and also for
    multiprecision objects.}
}
% \details{ _____ FIXME: Give exact formula /citation of Abramowitz & Stegun
% }                      and possibly also the DLMF formula + URL
\value{
  a \dQuote{numeric} (or complex or \code{"\link[Rmpfr:mpfr-class]{mpfr}"})
  vector of the same class and length as \code{x}.
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
(nus <- c(outer((0:3)/4, 1:5, `+`)))
stopifnot(
  all.equal(besselIs(1:10, 1), # our R code
            besselI (1:10, 1)) # internal C code w/ different algorithm
  ,
  sapply(nus, function(nu)
   all.equal(besselIs(1:10, nu, expon.scale=TRUE), # our R code
             BesselI (1:10, nu, expon.scale=TRUE)) # TOMS644 code
   )
  ,
  ## complex argument [gives warnings  'nterm=800' may be too small]
  sapply(nus, function(nu)
   all.equal(besselIs((1:10)*(1+1i), nu, expon.scale=TRUE), # our R code
             BesselI ((1:10)*(1+1i), nu, expon.scale=TRUE)) # TOMS644 code
   )
)

## Large 'nu' ...
x <- (0:20)/4
(bx <- besselI(x, nu=200))# base R's -- gives (mostly wrong) warnings
if(require("Rmpfr")) { ## Use high precision (notably large exponent range) numbers:
  Bx <- besselIs(mpfr(x, 64), nu=200)
  all.equal(Bx, bx, tol = 1e-15)# TRUE -- warning were mostly wrong; specifically:
  cbind(bx, Bx)
  signif(asNumeric(1 - (bx/Bx)[19:21]), 4) # only [19] had lost accuracy

  ## With*out* mpfr numbers -- using log -- is accurate (here)
  (lbx <- besselIs(     x,      nu=200, log=TRUE))
  lBx <-  besselIs(mpfr(x, 64), nu=200, log=TRUE)
  stopifnot(all.equal(asNumeric(log(Bx)), lbx, tol=1e-15),
	    all.equal(lBx, lbx, tol=4e-16))
} # Rmpfr
}
\keyword{math}
