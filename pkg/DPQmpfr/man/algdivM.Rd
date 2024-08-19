\name{algdivM}
\title{Compute log(gamma(b)/gamma(a+b)) Accurately, also via \pkg{Rmpfr}}
\alias{algdivM}
\description{
  Computes \deqn{\code{algdiv(a,b)} := \log \frac{\Gamma(b)}{\Gamma(a+b)} = \log
    \Gamma(b) - \log\Gamma(a+b) = \code{lgamma(b) - lgamma(a+b)}}{%
    \code{algdiv(a,b)} := log (\Gamma(b)/\Gamma(a+b)) =
    log \Gamma(b) - log\Gamma(a+b) = \code{lgamma(b) - lgamma(a+b)}}
  in a numerically stable way.

  The name \sQuote{algdiv} is from the auxiliary function in \R's (TOMS 708) implementation of
  \code{\link{pbeta}()}.
  As package \CRANpkg{DPQ} provides \R's Mathlib (double precision) as \R
  function \code{algdiv()}, we append \sQuote{M} to show the reliance on
  the \CRANpkg{Rmpfr} package.
}
\usage{% --> ../R/special-fun.R
algdivM(a, b, usePr = NULL)
}
\arguments{
  \item{a, b}{numeric or numeric-alike vectors (recycled to the same length
    if needed), typically inheriting from \code{class} \code{"mpfr"}.}
  \item{usePr}{positive integer specifying the precision in \bold{bit}s, or
    \code{\link{NULL}} when a smart default will be used.}
}
\details{
  Note that this is also useful to compute the Beta function
  \deqn{B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}.}{B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b).}
  Clearly,
  \deqn{\log B(a,b) = \log\Gamma(a) + \mathrm{algdiv(a,b)}
                    = \log\Gamma(a) - \mathrm{logQab}(a,b).}{%
         log B(a,b) = log\Gamma(a) + algdiv(a,b) = log\Gamma(a) - logQab(a,b).}

  In our \file{../tests/qbeta-dist.R} file, we look into computing
  \eqn{\log(p B(p,q))}{log(p * Beta(p,q))} accurately for
  \eqn{p \ll q}{p << q} .

  We are proposing a nice solution there. %% __ FIXME :
  \cr How is this related to \code{algdiv()} ?

  Additionally, we have defined
  \deqn{Qab = Q_{a,b} := \frac{\Gamma(a+b),\Gamma(b)},}{%
        Qab           :=       \Gamma(a+b)/\Gamma(b),}
  such that \eqn{\code{logQab(a,b)} := \log Qab(a,b)}{%
                 \code{logQab(a,b)} :=  log(Qab(a,b))}
  fulfills simply
  \deqn{\code{logQab(a,b)} = - \code{algdiv(a,b)}}
  see \code{\link[DPQ]{logQab_asy}} from package \CRANpkg{DPQ}.
}
\value{
  a numeric vector of length \code{max(length(a), length(b))} (if neither
  is of length 0, in which case the result has length 0 as well).
}
\references{
  Didonato, A. and Morris, A., Jr, (1992)
  Algorithm 708: Significant digit computation of the incomplete beta
  function ratios,
  \emph{ACM Transactions on Mathematical Software} \bold{18}, 360--373. % =: TOMS 708
}
\author{Martin Maechler (for the \pkg{Rmpfr} version).}
\seealso{
  \code{\link{gamma}}, \code{\link{beta}};
  the (double precision) version \code{\link[DPQ]{algdiv}()} in \CRANpkg{DPQ},
  and also in \pkg{DPQ}, the asymptotic approximation
  \code{\link[DPQ]{logQab_asy}()}.% <--> ~/R/Pkgs/DPQ/man/lbeta.Rd
}
\examples{
Qab <- algdivM(2:3, 8:14)
cbind(a = 2:3, b = 8:14, Qab) # recycling with a warning

## algdivM()  and my  logQab_asy()  give *very* similar results for largish b:
(lQab <- DPQ::logQab_asy(3, 100))
all.equal( - algdivM(3, 100), lQab, tolerance=0) # 1.283e-16 !!
## relative error
1 + lQab/ algdivM(3, 1e10) # 0 (64b F 30 Linux; 2019-08-15)

## in-and outside of "certified" argument range {b >= 8}:
a. <- c(1:3, 4*(1:8))/32
b. <- seq(1/4, 20, by=1/4)
ad <- t(outer(a., b., algdivM))
## direct computation:
f.algdiv0 <- function(a,b) lgamma(b) - lgamma(a+b)
f.algdiv1 <- function(a,b) lgamma(b) - lgamma(a+b)
ad.d <- t(outer(a., b., f.algdiv0))

matplot (b., ad.d, type = "o", cex=3/4,
         main = quote(log(Gamma(b)/Gamma(a+b)) ~"  vs.  algdivM(a,b)"))
mtext(paste0("a[1:",length(a.),"] = ",
        paste0(paste(head(paste0(formatC(a.*32), "/32")), collapse=", "), ", .., 1")))
matlines(b., ad,   type = "l", lwd=4, lty=1, col=adjustcolor(1:6, 1/2))
abline(v=1, lty=3, col="midnightblue")
# The larger 'b', the more accurate the direct formula wrt algdivM()
all.equal(ad[b. >= 1,], ad.d[b. >= 1,]       )# 1.5e-5
all.equal(ad[b. >= 2,], ad.d[b. >= 2,], tol=0)# 3.9e-9
all.equal(ad[b. >= 4,], ad.d[b. >= 4,], tol=0)# 4.6e-13
all.equal(ad[b. >= 6,], ad.d[b. >= 6,], tol=0)# 3.0e-15
all.equal(ad[b. >= 8,], ad.d[b. >= 8,], tol=0)# 2.5e-15 (not much better)
%% eventually use sfsmisc::relErr[V]() for visualization *and* checking
}
\keyword{math}
