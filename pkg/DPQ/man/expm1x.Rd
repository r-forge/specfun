\name{expm1x}
\title{Accurate exp(x) - 1 - x   (for smallish |x|)}
\alias{expm1x}
\alias{expm1xTser}
\description{
  Compute \eqn{e^x - 1 - x =} \code{exp(x) - 1 - x} accurately, notably for small \eqn{|x|}.

  The last two entries in \code{cutx[]} denote boundaries where
  \code{expm1x(x)} uses direct formulas.  For \code{nC <- length(cutx)},
  \code{exp(x) - 1 - x} is used for \code{abs(x) >= cutx[nC]},  and when \code{abs(x) < cutx[nC]}
  \code{expm1(x) - x}   is used for \code{abs(x) >= cutx[nC-1]}.
}
\usage{
expm1x(x, cutx = c( 4.4e-8, 0.1, 0.385, 1.1, 2),
             k = c(2,      9,  12,    17))

expm1xTser(x, k)
}
\arguments{
  \item{x}{numeric-alike vector; goal is to work for
    \code{\link[Rmpfr]{mpfr}}-numbers too.}
  \item{cutx}{increasing positive numeric vector of cut points defining
    intervals in which the computations will differ.}
  \item{k}{for \describe{
      \item{\code{exp1mx()}:}{increasing vector of integers with
	\code{\link{length}(k) == length(cutx) + 2}, denoting the order of
	Taylor polynomial approximation by \code{expm1xTser(.,k)} to \code{expm1x(.)}.}
      \item{\code{exp1mxTser()}:}{an integer \eqn{\ge 1}{>= 1}, where the
	Taylor polynomial approximation has degree \eqn{k + 1}.}
  }}
}
%% \details{
%% }
\value{
  a vector like \code{x} containing (approximations to) \eqn{e^x - x - 1}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{expm1}(x)} for computing \eqn{e^x - 1} is much more widely
  know, and part of the ISO C standards now.
}
\examples{%% --- more, notably tests: see ../tests/expm1x-tst.R

## a symmetric set of negative and positive
x <- unique(c(2^-seq(-3/8, 54, by = 1/8), seq(7/8, 3, by = 1/128)))
x <- x0 <- sort(c(-x, 0, x)) # negative *and* positive

## Mathematically,  expm1x() = exp(x) - 1 - x  >= 0  (and == 0 only at x=0):
em1x <- expm1x(x)
stopifnot(em1x >= 0, identical(x == 0, em1x == 0))

plot (x, em1x, type='b', log="y")
lines(x, expm1(x)-x, col = adjustcolor(2, 1/2), lwd = 3) ## should nicely cover ..
lines(x, exp(x)-1-x, col = adjustcolor(4, 1/4), lwd = 5) ## should nicely cover ..
cuts <- c(4.4e-8, 0.10, 0.385, 1.1, 2)[-1] # *not* drawing 4.4e-8
v <- c(-rev(cuts), 0, cuts); stopifnot(!is.unsorted(v))
abline(v = v, lty = 3, col=adjustcolor("gray20", 1/2))

stopifnot(diff(em1x[x <= 0]) <= 0)
stopifnot(diff(em1x[x >= 0]) >= 0)

## direct formula - may be really "bad" :
expm1x.0 <- function(x) exp(x) -1 - x
## less direct formula - improved (but still not universally ok):
expm1x.1 <- function(x) expm1(x)  - x

ax <- abs(x) # ==> show negative and positive x on top of each other
plot (ax, em1x, type='l', log="xy", xlab = "|x|  (for negative and positive x)")
lines(ax, expm1(x)-x, col = adjustcolor(2, 1/2), lwd = 3) ## see problem at very left
lines(ax, exp(x)-1-x, col = adjustcolor(4, 1/4), lwd = 5) ## see huge problems for |x| < ~10^{-7}
legend("topleft", c("expm1x(x)", "expm1(x) - x", "exp(x) - 1 - x"), bty="n",
       col = c(1,2,4), lwd = c(1,3,5))

## -------------------- Relative error of Taylor series approximations :
twoP <- seq(-0.75, 54, by = 1/8)
x <- 2^-twoP
x <- sort(c(-x,x)) # negative *and* positive
e1xAll <- cbind(expm1x.0 = expm1x.0(x),
                expm1x.1 = expm1x.1(x),
                vapply(1:15, \(k) expm1xTser(x, k=k), x))
colnames(e1xAll)[-(1:2)] <- paste0("k=",1:15)
head(e1xAll)
## TODO  plot !!
}
\keyword{math}
\keyword{arith}
