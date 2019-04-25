\name{Airy}
\Rdversion{1.1}
\alias{Airy}
\alias{AiryA}
\alias{AiryB}
\title{Airy Functions (and Their First Derivative)}
\description{
  Compute the Airy functions \eqn{Ai} or \eqn{Bi} or their first
  derivatives,
  \eqn{\frac{d}{dz} Ai(z)}{d/dz Ai(z)} and
  \eqn{\frac{d}{dz} Bi(z)}{d/dz Bi(z)}.

  The Airy functions are solutions of the differential equation
  \deqn{w'' = z w}
  for \eqn{w(z)}, and are related to each other and to the
  (modified) Bessel functions via (many identities, see
  \url{https://dlmf.nist.gov/9.6}), e.g., if
  \eqn{\zeta := \frac{2}{3} z \sqrt{z} = \frac{2}{3} z^{\frac{3}{2}}}{%
    zta := (2/3) z * sqrt(z) = 2/3 z^(3/2)},
  \deqn{Ai(z) = \pi^{-1}\sqrt{z/3}K_{1/3}(\zeta) =
            \frac{1}{3}\sqrt{z}\left(I_{-1/3}(\zeta) - I_{1/3}(\zeta)\right),}{%
    Ai(z) = 1/pi sqrt(z/3) K[1/3](zta) =
            1/3 sqrt(z) * (I[-1/3](zta) - I[1/3](zta)),}
 and
  \deqn{Bi(z) = \sqrt{z/3} \left(I_{-1/3}(\zeta) + I_{1/3}(\zeta)\right).}{%
    Bi(z) = sqrt(z/3) * (I[-1/3](zta) + I[1/3](zta)).}
}
\usage{
AiryA(z, deriv = 0, expon.scaled = FALSE, verbose = 0)
AiryB(z, deriv = 0, expon.scaled = FALSE, verbose = 0)
}
\arguments{
  \item{z}{complex or numeric vector.}
  \item{deriv}{order of derivative; must be 0 or 1.}
  \item{expon.scaled}{logical indicating if the result should be scaled
    by an exponential factor (typically to avoid under- or over-flow).}
  \item{verbose}{integer defaulting to 0, indicating the level of verbosity
    notably from C code.}
}
\details{
    By default, when \code{expon.scaled} is false, \code{AiryA()}
    computes the complex Airy function \eqn{Ai(z)} or its derivative
    \eqn{\frac{d}{dz} Ai(z)}{d/dz Ai(z)} on \code{deriv=0} or \code{deriv=1}
    respectively. \cr
    When \code{expon.scaled} is true, it returns
    \eqn{\exp(\zeta) Ai(z)}{exp(zta)*Ai(z)} or
    \eqn{\exp(\zeta) \frac{d}{dz} Ai(z)}{exp(zta)* d/dz Ai(z)},
    effectively removing the exponential decay in
    \eqn{-\pi/3 < \arg(z) < \pi/3}{-pi/3 < arg(z) < pi/3} and
    the exponential growth in
    \eqn{\pi/3 < \left|\arg(z)\right| < \pi}{pi/3 < abs(arg(z)) < pi},
    where \eqn{\zeta= \frac{2}{3} z \sqrt{z}}{zta=(2/3)*z*sqrt(z)}, and
    \eqn{\arg(z) = }{arg(z) = }\code{\link{Arg}(z)}.

    While the Airy functions \eqn{Ai(z)} and \eqn{d/dz Ai(z)} are
    analytic in the whole \eqn{z} plane, the corresponding scaled
    functions (for \code{expon.scaled=TRUE}) have a cut along the
    negative real axis.
% FIXME: Do we know an analogous fact for Bi() .. ??

    By default, when \code{expon.scaled} is false, \code{AiryB()}
    computes the complex Airy function \eqn{Bi(z)} or its derivative
    \eqn{\frac{d}{dz} Bi(z)}{d/dz Bi(z)} on \code{deriv=0} or \code{deriv=1}
    respectively. \cr
    When \code{expon.scaled} is true, it returns
    \eqn{exp(-\left|\Re(\zeta)\right|) Bi(z)}{exp(-abs(Re(zta)))*Bi(z)} or
    \eqn{exp(-\left|\Re(\zeta)\right|)\frac{d}{dz}Bi(z)}{%
         exp(-abs(Re(zta)))* dBi(z)/dz},
    to remove the exponential behavior in both the left and right half
    planes where, as above,
    \eqn{\zeta= \frac{2}{3}\cdot z \sqrt{z}}{zta=(2/3)*z*sqrt(z)}.
}
\value{
  a complex or numeric vector of the same length (and class) as \code{z}.
}
\references{see \code{\link{BesselJ}}; notably for many results the

    Digital Library of Mathematical Functions (DLMF), Chapter 9
    \emph{Airy and Related Functions} at \url{https://dlmf.nist.gov/9}.
}
\author{
  Donald E. Amos, Sandia National Laboratories, wrote the original
  fortran code.
  Martin Maechler did the \R interface.
}
\seealso{
  \code{\link{BesselI}} etc; the Hankel functions \code{\link{Hankel}}.

  The CRAN package \CRANpkg{Rmpfr} has \code{\link[Rmpfr]{Ai}(x)} for
  arbitrary precise "mpfr"-numbers \code{x}.
}
\examples{
## The AiryA() := Ai() function -------------

curve(AiryA, -20, 100, n=1001)
curve(AiryA,  -1, 100, n=1011, log="y") -> Aix
curve(AiryA(x, expon.scaled=TRUE), -1, 50, n=1001)
## Numerically "proving" the 1st identity above :
z <- Aix$x; i <- z > 0; head(z <- z[i <- z > 0])
Aix <- Aix$y[i]; zeta <- 2/3*z*sqrt(z)
stopifnot(all.equal(Aix, 1/pi * sqrt(z/3)* BesselK(zeta, nu = 1/3),
                    tol = 4e-15)) # 64b Lnx: 7.9e-16;  32b Win: 1.8e-15

\donttest{## This gives many warnings (248 on nb-mm4, F24) about lost accuracy, but on Windows takes ~ 4 sec:
curve(AiryA(x, expon.scaled=TRUE),  1, 10000, n=1001, log="xy")}

## The AiryB() := Bi() function -------------

curve(AiryB, -20, 2, n=1001); abline(h=0,v=0, col="gray",lty=2)
curve(AiryB, -1, 20, n=1001, log = "y") # exponential growth (x > 0)

curve(AiryB(x,expon.scaled=TRUE), -1, 20,    n=1001)
curve(AiryB(x,expon.scaled=TRUE),  1, 10000, n=1001, log="x")
}
\keyword{math}
