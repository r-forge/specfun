\name{stirlerrM}
\title{Stirling Formula Approximation Error}
\alias{stirlerrM}
\alias{stirlerrSer}
\description{
  Compute the \code{\link{log}()} of the error of Stirling's formula for \eqn{n!}.
  Used in certain accurate approximations of (negative) binomial and Poisson probabilities.

  \code{stirlerrM()} currently simply uses the direct mathematical formula,
  based on \code{\link{lgamma}()}, adapted for use with \code{\link[Rmpfr]{mpfr}}-numbers.
}
\usage{
stirlerrM(n, minPrec = 128L)
stirlerrSer(n, k)
}
\arguments{
  \item{n}{numeric or \dQuote{numeric-alike} vector, typically
    \dQuote{large} positive integer or half integer valued, here typically an
    \code{"mpfr"}-number vector.}
  \item{k}{integer \emph{scalar}, now in \code{1:11}.}
  \item{minPrec}{minimal precision (in bits) to be used when coercing
    number-alikes, say, biginteger (\code{\link[gmp]{bigz}}) to \code{"mpfr"}.}
}
\details{
  Stirling's approximation to \eqn{n!} has been
  \deqn{n! \approx \bigl(\frac{n}{e}\bigr)^n \sqrt{2\pi n},}{%
    n! ~= (n/e)^n * sqrt(2*pi*n),}
  where by definition the error is the difference of the left and right
  hand side of this formula, in  \eqn{\log}-scale,
  \deqn{\delta(n) =  \log\Gamma(n + 1) - n \log(n) + n - \log(2 \pi n)/2.}{%
         delta(n) =   log\Gamma(n + 1) - n* log(n) + n -  log(2*pi*n)/2.}

  See the vignette \emph{log1pmx, bd0, stirlerr, ...} from package
  \CRANpkg{DPQ}, where the series expansion of \eqn{\delta(n)} is used with
  11 terms, starting with
  \deqn{\delta(n) = \frac 1{12 n} - \frac 1{360 n^3} + \frac 1{1260 n^5}
                    \pm O(n^{-7}).}{%
	 delta(n) = 1/(12 n) - 1/(360 n^3) + 1/(1260 n^5) +/- O(n^{-7}).}
}
\note{
  In principle, the direct formula should be replaced by a few terms of the
  series in powers of \eqn{1/n} for large \code{n}, but we assume using
  high enough precision for \code{n} should be sufficient and \dQuote{easier}.
}
\value{
  a numeric or other \dQuote{numeric-alike} class vector, e.g.,
  \code{\link[Rmpfr]{mpfr}}, of the same length as \code{n}.
}
\references{
  Catherine Loader, see \code{\link{dbinom}};

  %% source: ~/R/Pkgs/DPQ/vignettes/log1pmx-etc.Rnw
  Martin Maechler (2021)
  log1pmx(), bd0(), stirlerr() -- Computing Poisson, Binomial, Gamma Probabilities in \R.
  \url{https://CRAN.R-project.org/package=DPQ/vignettes/log1pmx-etc.pdf}
}
\author{Martin Maechler}
\seealso{
 \code{\link{dbinom}}, \code{\link[DPQ]{stirlerr}()} in package
  \CRANpkg{DPQ} which is a pure \R version \R's mathlib-internal C function.
}
\examples{
### ----------------  Regular R  double precision -------------------------------

n <- n. <- c(1:10, 15, 20, 30, 50*(1:6), 100*(4:9), 10^(3:12))
(stE <- stirlerrM(n)) # direct formula is *not* good when n is large:
require(graphics)
plot(stirlerrM(n) ~ n, log = "x", type = "b", xaxt="n")
sfsmisc::eaxis(1, sub10=3)
for(k in 1:8) lines(n, stirlerrSer(n, k), col = k+1)
legend("top", c("stirlerrM(n)", paste0("stirlerrSer(n, k=", 1:8, ")")),
       pch=c(1,rep(NA,8)), col=1:(8+1), lty=1, bty="n")
## for larger n, current values are even *negative* ==> dbl prec *not* sufficient

## y in log-scale [same conclusion]
plot (stirlerrM(n) ~ n, log = "xy", type = "b", ylim = c(1e-13, 0.08))
for(k in 1:8) lines(n, stirlerrSer(n, k), col = k+1)
legend("topright", c("stirlerrM(n)", paste0("stirlerrSer(n, k=", 1:8, ")")),
       pch=c(1,rep(NA,8)), col=1:(8+1), lty=1, ncol=2, bty="n")

## the numbers:
options(digits=4, width=111)

stEmat. <- cbind(sM = stirlerrM(n),
                sapply(setNames(1:8, paste0("k=",1:8)),
                       function(k) stirlerrSer(n=n, k=k)))
stEmat.

## for printing n=<nice>:
N <- Rmpfr::asNumeric
dfm <- function(n, mm) data.frame(n=formatC(N(n)), N(mm), check.names=FALSE)

## relative differences:
dfm(n, stEmat.[,-1]/stEmat.[,1] - 1)
    # => stirlerrM() {with dbl prec} deteriorates after ~ n = 200--500
dfm(n, stEmat.[,-(1+8)]/stEmat.[,1+8] - 1)


### ----------------  MPFR High Accuracy -------------------------------

stopifnot(require(gmp),
          require(Rmpfr))
n <- as.bigz(n.)
## now repeat everything .. from above ... FIXME shows bugs !
## fully accurate using big rational arithmetic
class(stEserQ <- sapply(setNames(1:8, paste0("k=",1:8)),
                        function(k) stirlerrSer(n=n, k=k))) # list ..
stopifnot(sapply(stEserQ, class) == "bigq") # of exact big rationals
str(stEsQM  <- lapply(stEserQ, as, Class="mpfr"))# list of 8;  each prec. 128..702
    stEsQM. <- lapply(stEserQ, .bigq2mpfr, precB = 512) # constant higher precision
stEsQMm <- sapply(stEserQ, asNumeric) # a matrix

stEM   <- stirlerrM(mpfr(n, 128)) # now ok (loss of precision, but still ~ 10 digits correct)
stEM4k <- stirlerrM(mpfr(n, 4096))# assume "perfect"
## ==> what's the accuracy of the 128-bit 'stEM'?
N <- asNumeric # short
dfm(n, stEM/stEM4k - 1)
## 29 1e+06  4.470e-25
## 30 1e+07 -7.405e-23
## 31 1e+08 -4.661e-21
## 32 1e+09 -7.693e-20
## 33 1e+10  3.452e-17  (still ok)
## 34 1e+11 -3.472e-15  << now start losing
## 35 1e+12 -3.138e-13  <<<<
## same conclusion via  number of correct (decimal) digits:
dfm(n, log10(abs(stEM/stEM4k - 1)))

plot(N(-log10(abs(stEM/stEM4k - 1))) ~ N(n), type="o", log="x",
     xlab = quote(n), main = "#{correct digits} of 128-bit stirlerrM(n)")
ubits <- c(128, 52) # above 128-bit and double precision
abline(h = ubits* log10(2), lty=2)
text(1, ubits* log10(2), paste0(ubits,"-bit"), adj=c(0,0))

stopifnot(identical(stirlerrM(n), stEM)) # for bigz & bigq, we default to precBits = 128
all.equal(roundMpfr(stEM4k, 64),
          stirlerrSer (n., 8)) # 0.00212 .. because of 1st few n.  ==> drop these
all.equal(roundMpfr(stEM4k,64)[n. >= 3], stirlerrSer (n.[n. >= 3], 8)) # 6.238e-8

plot(asNumeric(abs(stirlerrSer(n., 8) - stEM4k)) ~ n.,
     log="xy", type="b", main="absolute error of stirlerrSer(n, 8)  & (n, 5)")
abline(h = 2^-52, lty=2); text(1, 2^-52, "52-bits", adj=c(1,-1)/8)
lines(asNumeric(abs(stirlerrSer(n., 5) - stEM4k)) ~ n., col=2)

plot(asNumeric(stirlerrM(n)) ~ n., log = "x", type = "b")
for(k in 1:8) lines(n, stirlerrSer(n, k), col = k+1)
legend("top", c("stirlerrM(n)", paste0("stirlerrSer(n, k=", 1:8, ")")),
       pch=c(1,rep(NA,8)), col=1:(8+1), lty=1, bty="n")

## y in log-scale
plot(asNumeric(stirlerrM(n)) ~ n., log = "xy", type = "b", ylim = c(1e-13, 0.08))
for(k in 1:8) lines(n, stirlerrSer(n, k), col = k+1)
legend("top", c("stirlerrM(n)", paste0("stirlerrSer(n, k=", 1:8, ")")),
       pch=c(1,rep(NA,8)), col=1:(8+1), lty=1, bty="n")
## all "looks" perfect (so we could skip this)

## the numbers ...
## %% FIXME a list instead of mpfrMatrix ... FIXME _____________
## FIXME ... asNumeric() needed or  as(*, "mpfr") or ...
ks <- 1:8 ## k <= 5 === FIXME --- use DPQ's  version !!
stirlS.l <- lapply(setNames(ks, paste0("k=",ks)),
                   function(k) stirlerrSer(n=n, k=k))
## ==> an mpfrMatrix of dim  35 x 5 :
mss <- do.call(cbind, lapply(stirlS.l, mpfr, precBits=256))
stEmat <- cbind(sM = stEM4k, mss)
signif(asNumeric(stEmat), 6) # so it prints nicely
## print *relative errors* nicely :
## simple double precision version of direct formula (cancellation for n >> 1 !):
stE <- stirlerrM(n.)
dfm(n , cbind(stEmat[,-1], dbl=stE)/stEM4k - 1)
## relative differences:
dfm(n,  stEmat[,-1]    / stEmat[,1] - 1)
dfm(n., stEmat[,-(1+8)]/ stEmat[,1+8] - 1)
}
\keyword{math}
\keyword{arith}
