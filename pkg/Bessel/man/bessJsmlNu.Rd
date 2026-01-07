\name{bessJsmlNu}
\alias{bessJsmlNu}
\title{Bessel \eqn{J_\nu(x)} for Very Small Order \eqn{\nu}}
\description{
  Use a few terms Taylor approximation for Bessel's \eqn{J_\nu(x)} for
  fixed \eqn{x} and \eqn{\nu \to 0}{\nu --> 0}.
}
\usage{
bessJsmlNu(x, nu, nTrms, m.max)
}
\arguments{
  \item{x, nu}{numeric Bessel J() arguments, see e.g., \code{\link{besselJ}}.}
  %% FIXME/TODO: complex arguments should work too!
  \item{nTrms}{integer, one of 0, 1, or 2, denoting the number of Taylor
    series terms to be used,  \code{nTrms = 0} simply returning
    \code{\link{BesselJ}(x, 0)}.}
  \item{m.max}{positive integer specifying the Taylor series order used.}
}
\details{
  A. & S.  (9.1.10) p. 360  (Ascending Series)
           (9.1.64) p. 362  (derivative wrt order)

  ............ show derivation? ...........
}
\value{
  The function explicitly calls \code{\link{Vectorize}(<fn>, c("x", "nu"))},
  hence a numeric vector conforming to (e.g., ) \code{x + nu}.
}
\references{
  Abramowitz, M., and Stegun, I. A. (1964--1972).
  \emph{Handbook of mathematical functions}
  (NBS AMS series 55, U.S. Dept. of Commerce).
  \url{https://personal.math.ubc.ca/~cbm/aands/page_360.htm}
}
\author{Martin Maechler}
%% \note{
%%   (note)
%% }
\seealso{
  This package's \code{\link{BesselJ}}, base \R's \code{\link{besselJ}}
}
\examples{
(nuSml <- 2^-c(seq(30, 53, by=1/2), 75, 100, 300, 800, 1022, 1074, Inf))
## 9.3e-10 6.6e-10 ... 1.1e-16 .. 2.22e-308 4.94e-324 0.0
## bug in R (at least up to Jan. 1 2026) :
options(digits = 14) # show more digits
Jsml <- sapply(nuSml, function(nu) besselJ(pi/2, nu))
tail(cbind(nuSml, Jsml), 20) # complete "divergence" for  nu <~ 10^-{16}
Jsm.nT.1 <- bessJsmlNu(pi/2, nuSml, nTrms = 1, m.max = 30)
Jsm.nT.2 <- bessJsmlNu(pi/2, nuSml, nTrms = 2, m.max = 30)
cbind(nuSml, Jsml, Jsm.nT.1, Jsm.nT.2)
table(Jsm.nT.2 == Jsm.nT.1) # all TRUE  (really ???)
          all.equal(Jsm.nT.2, Jsm.nT.1, tolerance = 0) # show
stopifnot(all.equal(Jsm.nT.2, Jsm.nT.1, tolerance = 4e-16))

## second example;  nu  *not* very small
x. <- 2
nu2 <- 2^-seq(1, 30, by = 1/8)
Jsm2.1 <- bessJsmlNu(x., nu2, nTrms = 1, m.max = 40)
Jsm2.2 <- bessJsmlNu(x., nu2, nTrms = 2, m.max = 40)
table(Jsm2.2 == Jsm2.1) # now they *do* differ
bessJ2 <- besselJ(x., nu2)
BessJ2 <- sapply(nu2, function(nu) BesselJ(x., nu))
## J() function values: ---------------------
matplot(nu2, cbind(bessJ2, BessJ2, Jsm2.1, Jsm2.2), log = "x", type = "l")
sum(i <- nu2 < 1e-4)
stopifnot(all.equal(BessJ2[i], Jsm2.2[i], tolerance = 4e-14),
          all.equal(BessJ2[i], Jsm2.1[i], tolerance = 4e-9))

##  |error|  wrt BesselJ() ---------------
matplot(nu2, abs(cbind(bessJ2, Jsm2.1, Jsm2.2) - BessJ2), log = "xy",
        type = "l", ylim = c(5e-17, 1e-10),
        xlab = quote("nu" == nu), xaxt = "n", yaxt = "n",
	main = substitute(list(abs("<bess J>"(x, nu) - BesselJ(x, nu)), x == XX),
                          list(XX = formatC(x.))))
eaxis <- if(!requireNamespace("sfsmisc")) axis else sfsmisc::eaxis
eaxis(1, sub10 = c(-2,0)); eaxis(2)
legend("topleft", lty = 1:5, lwd = 2, col = 1:6, bty = "n",
       legend = paste0(c("besselJ(", paste0("bessJsmlNu(*, nTrms = ",1:2)), ")"))
}
\keyword{math}
