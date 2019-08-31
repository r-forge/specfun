\name{lssum}
\alias{lssum}
\title{Compute Logarithm of a Sum with Signed Large Summands}
\usage{
lssum(lxabs, signs, l.off = max(lxabs), strict = TRUE)
}
\arguments{
  \item{lxabs}{n-vector of values \eqn{\log(|x_1|), \ldots, \log(|x_n|)}.}
  \item{signs}{corresponding signs \eqn{sign(x_1), \ldots, sign(x_n)}.}
  \item{l.off}{the offset to substract and re-add; ideally in the order of \code{\link{max}(.)}.}
  \item{strict}{\code{\link{logical}} indicating if the function should stop on some negative sums.}
}
\value{
  \deqn{
    log(x_1 + .. + x_n) =
    = log(sum(x)) = log(sum(sign(x)*|x|)) =
    = log(sum(sign(x)*exp(log(|x|)))) =
    = log(exp(log(x0))*sum(signs*exp(log(|x|)-log(x0)))) =
    = log(x0) + log(sum(signs* exp(log(|x|)-log(x0)))) =
    = l.off + log(sum(signs* exp(lxabs - l.off )))
  }
}
\description{
  Properly compute \eqn{\log(x_1 + \ldots + x_n)}{log(x1 + .. + xn)}
  for given log absolute values \code{lxabs = }
  \eqn{log(|x_1|),.., log(|x_n|)}{log(|x_1|),.., log(|x_n|)}
  and corresponding signs \code{signs = } \eqn{sign(x_1),.., sign(x_n)}.  Here,
  \eqn{x_i}{xi} is of arbitrary sign.

  Notably this works in many cases where the direct sum would have summands
  that had overflown to \code{+Inf} or underflown to \code{-Inf}.

  This is a (simpler, vector-only) version of \code{copula:::lssum()} (CRAN
  package \CRANpkg{copula}). %  ~/R/Pkgs/copula/R/special-func.R

  Note that the \emph{precision} is often not the problem for the direct
  summation, as \R's \code{\link{sum}()} internally uses
  \code{"long double"} precision on most platforms.
}
\seealso{
  \code{\link{lsum}()} which computes an exponential sum in log scale
  with\emph{out} signs.
}
\author{
  Marius Hofert and Martin Maechler (for package \CRANpkg{copula}).
}
\examples{
rSamp <- function(n, lmean, lsd = 1/4, roundN = 16) {
  lax <- sort((1+1e-14*rnorm(n))*round(roundN*rnorm(n, m = lmean, sd = lsd))/roundN)
  sx <- rep_len(c(-1,1), n)
  list(lax=lax, sx=sx, x = sx*exp(lax))
}

set.seed(101)
L1 <- rSamp(1000, lmean = 700) # here, lssum() is not needed (no under-/overflow)
summary(as.data.frame(L1))
ax <- exp(lax <- L1$lax)
hist(lax); rug(lax)
hist( ax); rug( ax)
sx <- L1$sx
table(sx)
(lsSimple <- log(sum(L1$x)))           # 700.0373
(lsS <- lssum(lxabs = lax, signs = sx))# ditto
lsS - lsSimple # even exactly zero (in 64b Fedora 30 Linux which has nice 'long double')
stopifnot(all.equal(700.037327351478, lsS, tol=1e-14), all.equal(lsS, lsSimple))

L2 <- within(L1, { lax <- lax + 10; x <- sx*exp(lax) }) ; summary(L2$x) # some -Inf, +Inf
(lsSimpl2 <- log(sum(L2$ x)))                    # NaN
(lsS2 <- lssum(lxabs = L2$ lax, signs = L2$ sx)) # 710.0373
stopifnot(all.equal(lsS2, lsS + 10, tol = 1e-14))
}
\keyword{distribution}
