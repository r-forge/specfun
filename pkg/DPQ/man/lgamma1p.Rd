\name{lgamma1p}
\title{Accurate \code{log(gamma(a+1))}}
\alias{lgamma1p}
\alias{lgamma1p.}
\alias{lgamma1pC}
\alias{lgamma1p_series}
\description{
  Compute
  \deqn{l\Gamma_1(a) := \log\Gamma(a+1) = \log(a\cdot \Gamma(a)) = \log a + \log \Gamma(a),}{%
        lG[1](a) := log(Gamma(a+1)) = log(a * Gamma(a)) = log(a) + log(Gamma(a)),}
  which is \dQuote{in principle} the same as
  \code{\link{log}(\link{gamma}(a+1))} or \code{\link{lgamma}(a+1)},
  accurately also for (very) small \eqn{a} \eqn{(0 < a < 0.5)}.
}
\usage{
lgamma1p (a, tol_logcf = 1e-14, f.tol = 1, ...)	% ../R/beta-fns.R
lgamma1p.(a, cutoff.a = 1e-6, k = 3)		% ../R/qchisqAppr.R
lgamma1p_series(x, k)               		% ../R/beta-fns.R
lgamma1pC(x)                        		% ../R/utils.R
}
\arguments{
  \item{a, x}{a numeric vector.}
  \item{tol_logcf}{for \code{lgamma1p()}: a non-negative number passed to
    \code{\link{logcf}()} (and \code{\link{log1pmx}()} which calls \code{logcf()}).}
  \item{f.tol}{numeric (\bold{f}actor) used in
    \code{\link{log1pmx}(*, tol_logcf = f.tol * tol_logcf)}.}
  \item{...}{further optional arguments passed on to \code{\link{log1pmx}()}.}
  \item{cutoff.a}{for \code{lgamma1p.()}: a positive number indicating
    the cutoff to switch from  ...}% FIXME
  \item{k}{an integer, the number of terms in the series expansion used
    internally; currently for \describe{
      \item{\code{lgamma1p.()}: }{\eqn{k \le 3}{k <= 3}}
      \item{\code{lgamma1p_series()}:}{\eqn{k \le 15}{k <= 15}}
  }}
}
\details{
  \code{lgamma1p()} is an \R translation of the function (in Fortran) in
  Didonato and Morris (1992) which uses a 40-degree polynomial approximation.

  \code{lgamma1p.(u)} for small \eqn{|u|} uses up to 4 terms of
  \deqn{\Gamma(1+u) = 1 + u*(-\gamma_E + a_0 u + a_1 u^2 + a_2 u^3) + O(u^5),}
  where \eqn{a_0 := (\psi'(1) + \psi(1)^2)/2 = (\pi^2/6 + \gamma_E^2)/2},
  and \eqn{a_1} und \eqn{a_2} are similarly determined. % --> ../R/qchisqAppr.R
  Then \code{log1p(.)} of the \eqn{\Gamma(1+u) - 1} approximation above is used.

  \code{lgamma1p_series(x, k)} is a Taylor series approximation of order
  \code{k}, directly of \eqn{l\Gamma_1(a) := \log \Gamma(a+1)} (mostly via
  \command{Maple}), which starts as
  \eqn{-\gamma_E x + \pi^2 x^2/ 12 + \dots}{gammaE*x + pi^2 * x^2/ 12 + ...},
  where \eqn{\gamma_E}{gammaE} is Euler's constant 0.5772156649. % ...% TODO MORE  write!...

  \code{lgamma1pC()} is an interface to \R's C API (\file{Mathlib} / \file{Rmath.h})
  function \code{lgamma1p()}.
}
\value{
  a numeric vector with the same attributes as \code{a}.
}
\references{
  Didonato, A. and Morris, A., Jr, (1992)
  Algorithm 708: Significant digit computation of the incomplete beta function ratios.
  \emph{ACM Transactions on Mathematical Software}, \bold{18}, 360--373;
  see also \code{\link{pbeta}}.
}
\author{Morten Welinder (C code of Jan 2005, see R's bug issue
  \PR{7307}) for \code{lgamma1p()}.

  Martin Maechler, notably for \code{lgamma1p_series()} which works
  with package \CRANpkg{Rmpfr} but otherwise may be \emph{much} less
  accurate than Morten's 40 term series!
}
%% \note{
%% }
\seealso{
  \code{\link{log1pmx}}, \code{\link{log1p}}, \code{\link{pbeta}}.
}
\examples{
curve(lgamma1p, -1.25, 5, n=1001, col=2, lwd=2)
abline(h=0, v=-1:0, lty=c(2,3,2), lwd=c(1, 1/2,1))
for(k in 1:15)
  curve(lgamma1p_series(x, k=k), add=TRUE, col=adjustcolor(paste0("gray",25+k*4), 2/3), lty = 3)

curve(lgamma1p, -0.25, 1.25, n=1001, col=2, lwd=2)
abline(h=0, v=0, lty=2)
for(k in 1:15)
  curve(lgamma1p_series(x, k=k), add=TRUE, col=adjustcolor("gray20", 2/3), lty = 3)

curve(-log(x*gamma(x)), 1e-30, .8, log="xy", col="gray50", lwd = 3,
      axes = FALSE, ylim = c(1e-30,1)) # underflows to zero at x ~= 1e-16
eaxGrid <- function(at.x = 10^(1-4*(0:8)), at.y = at.x) {
    sfsmisc::eaxis(1, sub10 = c(-2, 2), nintLog=16)
    sfsmisc::eaxis(2, sub10 = 2, nintLog=16)
    abline(h = at.y, v = at.x, col = "lightgray", lty = "dotted")
}
eaxGrid()
curve(-lgamma( 1+x), add=TRUE, col="red2", lwd=1/2)# underflows even earlier
curve(-lgamma1p (x), add=TRUE, col="blue") -> lgxy
curve(-lgamma1p.(x), add=TRUE, col=adjustcolor("forest green",1/4),
      lwd = 5, lty = 2)
for(k in 1:15)
  curve(-lgamma1p_series(x, k=k), add=TRUE, col=paste0("gray",80-k*4), lty = 3)
stopifnot(with(lgxy, all.equal(y, -lgamma1pC(x))))

if(requireNamespace("Rmpfr")) { # accuracy comparisons, originally from  ../tests/qgamma-ex.R
    x <- 2^(-(500:11)/8)
    x. <- Rmpfr::mpfr(x, 200)
    ## versions of lgamma1p(x) := lgamma(1+x)
    ## lgamma1p(x) = log gamma(x+1) = log (x * gamma(x)) = log(x) + lgamma(x)
    xct. <- log(x.  * gamma(x.)) # using  MPFR  arithmetic .. no overflow/underflow ...
    xc2. <- log(x.) + lgamma(x.) #  (ditto)

    AllEq <- function(target, current, ...)
        Rmpfr::all.equal(target, current, ...,
                         formatFUN = function(x, ...) Rmpfr::format(x, digits = 9))
    print(AllEq(xct., xc2., tol = 0)) # 2e-57
    rr <- vapply(1:15, function(k) lgamma1p_series(x, k=k), x)
    colnames(rr) <- paste0("k=",1:15)
    relEr <- Rmpfr::asNumeric(sfsmisc::relErrV(xct., rr))
    ## rel.error of direct simple computation:
    relE.D <- Rmpfr::asNumeric(sfsmisc::relErrV(xct., lgamma(1+x)))

    matplot(x, abs(relEr), log="xy", type="l", axes = FALSE,
            main = "|rel.Err(.)| for lgamma(1+x) =~= lgamma1p_series(x, k = 1:15)")
    eaxGrid()
    p2 <- -(53:52); twp <- 2^p2; labL <- lapply(p2, function(p) substitute(2^E, list(E=p)))
    abline(h = twp, lty=3)
    axis(4, at=twp, las=2, line=-1, labels=as.expression(labL), col=NA,col.ticks=NA)
    legend("topleft", paste("k =", 1:15), ncol=3, col=1:6, lty=1:5, bty="n")
    lines(x, abs(relE.D), col = adjustcolor(2, 2/3), lwd=2)
    legend("top", "lgamma(1+x)", col=2, lwd=2)

    ## zoom in:
    matplot(x, abs(relEr), log="xy", type="l", axes = FALSE,
            xlim = c(1e-5, 0.1), ylim = c(1e-17, 1e-10),
            main = "|rel.Err(.)| for lgamma(1+x) =~= lgamma1p_series(x, k = 1:15)")
    eaxGrid(10^(-5:1), 10^-(17:10))
    abline(h = twp, lty=3)
    axis(4, at=twp, las=2, line=-1, labels=as.expression(labL), col=NA,col.ticks=NA)
    legend("topleft", paste("k =", 1:15), ncol=3, col=1:6, lty=1:5, bty="n")
    lines(x, abs(relE.D), col = adjustcolor(2, 2/3), lwd=2)
    legend("right", "lgamma(1+x)", col=2, lwd=2)

} # Rmpfr only

}
\keyword{distribution}
\keyword{math}
