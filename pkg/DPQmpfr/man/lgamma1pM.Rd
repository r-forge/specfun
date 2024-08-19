\name{lgamma1pM}
\title{Compute  log( Gamma(x+1) ) Arbitrarily (MPFR) Accurately}
\alias{gamln1}% name used in TOMS 708 and pkg DPQ (has *both* lgamma1p() & gamln1()
\alias{lgamma1pM}
\description{
  Computes \eqn{\log \Gamma(x+1)} accurately notably when \eqn{|x| \ll 1}{|x| << 1}.
  %% --> ../R/special-fun.R  and R`s  src/nmath/toms708.c
  For \code{"mpfr"} numbers, the precision is increased intermediately such
  that \eqn{a+1} should not lose precision.

  \R's "own" double prec version is
  soon % then "now" / "since 2024-..." / "from DPQ version ....."
  available in package in \CRANpkg{DPQ},
  %_ once gamln1()'s DPQ is on CRAN
  %_ under the name \code{\link[DPQ]{gamln1}()} (from TOMS 708).
  under the name \code{gamln1()} (from TOMS 708).
  %% e.g. ~/R/Pkgs/DPQ/man/gamln1.Rd

   %% FIXME1: R-only implementation is in
   %%                 ~/R/Pkgs/DPQ/TODO_R_versions_gam1_etc.R
}
\usage{% --> ../R/special-fun.R
lgamma1pM(a, usePr = NULL, DPQmethod = c("lgamma1p", "algam1"))
}
\arguments{
  \item{a}{a numeric or numeric-alike vector, typically inheriting from
    \code{class} \code{"mpfr"}.}
  \item{usePr}{positive integer specifying the precision in \bold{bit}s, or
    \code{\link{NULL}} when a smart default will be used.}
  \item{DPQmethod}{a character string; must be the name of an
    \code{lgamma1p()}-alike function from package \CRANpkg{DPQ}.  It will
    be called in case of \code{is.numeric(a)} (and when \pkg{DPQ} is available).}
}
%% \details{

%% }
\value{
  a numeric-alike vector like \code{a}.
}
\references{
  TOMS 708, see \code{\link{pbeta}}
}
\author{Martin Maechler}
%% \note{

%% }
\seealso{
  \code{\link{lgamma}()} (and \code{gamma()} (same page)),
  and our \code{\link{algdivM}()}; further, package \CRANpkg{DPQ}'s
  \code{\link[DPQ]{lgamma1p}()} and
  (if already available) \code{gamln1()}. %   \code{\link[DPQ]{gamln1}()}.
  %% NOT YET \code{\link{gam1}}, % NOT YET __FIXME__ \code{\link{pbeta_ser}}
}
\examples{
## Package {DPQ}'s  lgamma1p():
lgamma1p <- DPQ::lgamma1p
lg1 <- function(u) lgamma(u+1) # the simple direct form
u <- seq(-.5, 1.5, by=1/16); set.seed(1); u <- sample(u) # permuted (to check logic)
%% u <- seq(-.5, 1.5, by=1/4) ## for testing
g11   <- vapply(u, lgamma1p, numeric(1))
lgamma1p. <- lgamma1p(u)
all.equal(lg1(u), g11, tolerance = 0) # see 3.148e-16
stopifnot(exprs = {
    all.equal(lg1(u), g11, tolerance = 2e-15)
    identical(g11, lgamma1p.)
})

## Comparison using Rmpfr; slightly extending the [-.5, 1.5] interval:
u <- seq(-0.525, 1.525, length.out = 2001)
lg1p  <- lgamma1pM(   u)
lg1pM <- lgamma1pM(Rmpfr::mpfr(u, 128))
asNumeric <- Rmpfr::asNumeric
relErrV   <- sfsmisc::relErrV
if(FALSE) { # DPQ "latest" version __FIXME__
lng1  <- DPQ::lngam1(u)
relE <- asNumeric(relErrV(lg1pM, cbind(lgamma1p = lg1p, lngam1 = lng1)))
} else {
relE <- asNumeric(relErrV(lg1pM, cbind(lgamma1p = lg1p)))#, lngam1 = lng1)))
}

## FIXME: lgamma1p() is *NOT* good around u =1. -- even though it should
##        and the R-only vs (not installed) *does* "work" (is accurate there) ?????
## --> ~/R/Pkgs/DPQ/TODO_R_versions_gam1_etc.R
if(FALSE) {
matplot(u, relE, type="l", ylim = c(-1,1) * 2.5e-15,
     main = expression("relative error of " ~~ lgamma1p(u) == log( Gamma(u+1) )))
} else {
plot(relE ~ u, type="l", ylim = c(-1,1) * 2.5e-15,
     main = expression("relative error of " ~~ lgamma1p(u) == log( Gamma(u+1) )))
}
grid(lty = 3); abline(v = c(-.5, 1.5), col = adjustcolor(4, 1/2), lty=2, lwd=2)

## what about the direct formula -- how bad is it really ?
relED <- asNumeric(relErrV(lg1pM, lg1(u)))
lines(relED ~ u, col = adjustcolor(2, 1/2), lwd = 2)
}
\keyword{math}
