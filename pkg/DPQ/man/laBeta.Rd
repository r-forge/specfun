\name{laBeta}
\alias{laBeta}
\title{Accurate Approximation of \code{log(a * beta(a,b))} for small \code{a}}
\description{
  Uses a Taylor series to compute \code{log(a * \link{beta}(a,b))}\eqn{ = % no \Beta, use B
    \log(a B(a,b)) = \log a + \log B(a,b) = }\code{log(a) + lbeta(a,b)},
  where both direct formulas suffer from cancellation for small \eqn{a},
  i.e., already for \eqn{a < \frac 1 2} (slightly)
}
\usage{
laBeta(a, b, nT)
}
\arguments{
  \item{a, b}{numeric(-alike) vectors of positive numbers.  Needs \eqn{a < 1}
    for convergence.}
  \item{nT}{an integer \eqn{\ge 1}{>= 1}, the number of Taylor series terms
    to use for the approximation.  Currently without a default; typically
    \code{nT = 17} is enough to get a more values than with the direct formula.}
}
%% \details{

%% }
\value{
  a numeric vector (of \dQuote{confirming length}) with an approximation to
  \code{log(a * \link{beta}(a,b))}.
}
\references{
  The author derived the Taylor series by integrating
  % see my comments in ../R/beta-fns.R  --- find my hand-written notes !!!
  Abramowitz & Stegun, p.259, (6.3.14) series for \eqn{\psi(1+z)} to get a
  series for \eqn{log \Gamma(z+1) = \log(z \Gamma(z)) = \log z + \log\Gamma(z)}.
  A. & St.: see reference e.g., in \code{\link{lbeta_asy}}.
}
\author{Martin Maechler, January 2022, then October 2025}
%% \note{
%% }
\seealso{
  \R's beta function and its logarithm \code{\link{beta}()} and \code{lbeta()}.
}
\examples{
b <- c(.25, .5, 1:8, 12, 20, 49)
laB0 <- function(a,b) log(a * beta(a,b))
laB1 <- function(a,b) log(a) + lbeta(a,b)
a <- 1/32
labMat <- cbind(laB0= laB0(a, b), laB1 = laB1(a, b),
                laBe05 = laBeta(a=a, b = b, nT=5),
                laBe10 = laBeta(a=a, b = b, nT=10),
                laBe15 = laBeta(a=a, b = b, nT=15),
                laBe17 = laBeta(a=a, b = b, nT=17))
cbind(b, labMat)
stopifnot(exprs = {
    all.equal(labMat[,1], labMat[,"laB1"],   tolerance = 4e-14) # see 1.777e-14
    all.equal(labMat[,2], labMat[,"laBe05"], tolerance = 1e-06) # see 6.387e-7
    all.equal(labMat[,2], labMat[,"laBe10"], tolerance = 2e-11) # 1.037e-11
    all.equal(labMat[,2], labMat[,"laBe15"], tolerance = 8e-15) # 2.52e-15
})
cbind(b, laB1 = laB1(2,b), laBe1 = laBeta(2,b, nT=1), laBe2 = laBeta(2,b, nT=2),
      laBe20 = laBeta(2,b, nT=20))
if(requireNamespace("Rmpfr")) withAutoprint({
    asNumeric <- Rmpfr::asNumeric
    mpfr <- Rmpfr::mpfr
    beta <- Rmpfr::beta
  (laBM <- asNumeric(log(a * beta(mpfr(a,128),b)))) # the "true" value
  stopifnot(identical(laBM, asNumeric(laB0(mpfr(a,128), b))))
  relErr <- sfsmisc::relErr
  relErr(laBM, log(a * beta(a,b))) # 1.79e-14
  cbind(apply(labMat, 2, relErr, target = laBM))
  ## laB0   1.795634e-14
  ## laB1   2.189151e-15
  ## laBe05 6.387111e-07
  ## laBe10 1.037178e-11
  ## laBe15 3.787459e-16
  ## laBe17 1.666482e-16

  a <- 0.001; aM <- mpfr(a, 128) ## -- smaller a -- losing 3 decimals w/ direct formula:
  relErr(asNumeric(laB0(aM, b)) -> laBM, laB0(a,b))   # 5.55e-13
  relErr(laBM,                           laB1(a,b))   # 1.49e-13 slightly better
  relErr(laBM, laBeta(a,b, nT =  5)) # 2.33e-14 (better!)
  relErr(laBM, laBeta(a,b, nT = 10)) # 1.60e-16
  relErr(laBM, laBeta(a,b, nT = 20)) #  ditto, i.e., not better
})
}
\keyword{math}
