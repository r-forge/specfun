\name{logcf}
\title{Continued Fraction Approximation of Log-Related Power Series}
\alias{logcf}
\alias{logcfR_vec}
\alias{logcfR}
\description{
  Compute a continued fraction approximation to the series (infinite sum)
  \deqn{\sum_{k=0}^\infty \frac{x^k}{i +k\cdot d} = \frac{1}{i} + \frac{x}{i+d} +
    \frac{x^2}{i+2*d} + \frac{x^3}{i+3*d} + \ldots}{%
    sum(k=0,...,Inf; x^k/(i+k*d)) = 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ...}

   Needed as auxiliary function in \code{\link{log1pmx}()} and \code{\link{lgamma1p}()}.
}
\usage{					%--> ../R/beta-fns.R <<<
logcf (x, i, d, eps, maxit = 10000L, trace = FALSE)
logcfR(x, i, d, eps, maxit = 10000L, trace = FALSE)
logcfR_vec(x, i, d, eps, maxit = 10000L, trace = FALSE)
}
\arguments{
  \item{x}{numeric vector of values less than 1.
    "mpfr"-numbers (of potentially high precision, package \CRANpkg{Rmpfr}) work in
    \code{logcfR*(x,*)}.}
  \item{i}{positive numeric}
  \item{d}{non-negative numeric}
  \item{eps}{positive number, the convergence tolerance.}
  \item{maxit}{a positive integer, the maximal number of iterations or
    terms in the truncated series used.}
  \item{trace}{logical (or non-negative integer in the future) indicating
    if (and how much) diagnostic output should be printed to the console
    during the computations.}
}
\details{
  \describe{
    \item{\code{logcfR()}:}{a pure \R version where the iterations happen
      vectorized in \code{x}, only for those components \code{x[i]} they
      have not yet converged.  This is particularly beneficial for
      not-very-short \code{"mpfr"} vectors \code{x}, and still conceptually
      equivalent to the \code{logcfR_vec()} version.}
    \item{\code{logcfR_vec()}:}{a pure \R version where each \code{x[i]} is
      treated separately, hence \dQuote{properly} vectorized, but slowly so.}
    \item{\code{logcf()}:}{only for \code{\link{numeric}} \code{x}, calls
      into (a clone of) \R's own (non-API currently) \code{logcf()} C
      Rmathlib function.}
  }% describe
}
\value{
  a numeric-alike vector with the same attributes as \code{x}.  For the
   \code{logcfR*()} versions, an \code{"mpfr"} vector if \code{x} is one.
}
\author{Martin Maechler, based on \R's \file{nmath/pgamma.c} implementation.}
\note{
 Rescaling is done by (namespace hidden) \dQuote{global}
 \code{scalefactor} which is \eqn{2^{256}}, represented exactly (in
   \code{\link{double}} precision).
}
\seealso{
  \code{\link{lgamma1p}}, \code{\link{log1pmx}}, and
  \code{\link{pbeta}}, whose prinicipal algorithm has evolved from TOMS 708.
}
\examples{% considerably more in ../tests/dnbinom-tst.R
x <- (-2:1)/2
logcf (x, 2,3, eps=1e-7, trace=TRUE) # shows iterations for each x[]
logcfR_vec(x, 2,3, eps=1e-7, trace=TRUE) # 1 line per x[]
logcfR_vec(x, 2,3, eps=1e-7, trace= 2  ) # shows iterations for each x[]

n <- 2049; x <- seq(-1,1, length.out = n)[-n] ; stopifnot(diff(x) == 1/1024)
plot(x, (lcf <- logcf(x, 2,3, eps=1e-12)), type="l", col=2)
lcR <- logcfR_vec (x, 2,3, eps=1e-12); all.equal(lcf, lcR , tol=0)
lcR.<- logcfR(x, 2,3, eps=1e-12); all.equal(lcf, lcR., tol=0)
all.equal(lcR, lcR., tol=0) # TRUE
all.equal(lcf, lcR., tol=0) # TRUE (x86_64, Lnx)
stopifnot(exprs = {
  all.equal(lcf, lcR., tol=1e-14)# seen 0 (!)
  all.equal(lcR, lcR., tol=5e-16)# seen 0 above
})

l32 <- curve(logcf(x, 3,2, eps=1e-7), -3, 1, n = 1000)
abline(h=0,v=1, lty=3, col="gray50")
##
plot(y~x, l32, log="y", type="l", main= "logcf(x, i, d)  in log-scale",
     ylim = c(.01, 10), col=2); abline(v=1, lty=3, col="gray50")
## other (i,d) than the (3,2) needed for log1pmx():
curve(logcfR(x,  5,  4,  eps = 1e-5), n=1000, add=TRUE, col = 4)
curve(logcfR(x, 20,  1,  eps = 1e-5), n=1000, add=TRUE, col = 5)
curve(logcfR(x, 1/4,1/2, eps = 1e-5), n=1000, add=TRUE, col = 6)
i.d <- cbind(c(2,3), c(5,4), c(20,1), c(1/4, 1/2))
legend("topleft", apply(i.d, 2, \(k) paste0("(i=",k[1],", d=", k[2],")")),
       lwd = 2, col = c(2,4:6), bty="n")
}
\keyword{math}
