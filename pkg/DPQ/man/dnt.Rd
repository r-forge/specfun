\name{dnt}
\title{Non-central t-Distribution Density - Algorithms and Approximations}
\alias{dntJKBf1}
\alias{dntJKBf}
\alias{.dntJKBch1}
\alias{.dntJKBch}
\description{
  \code{dntJKBf1} implements the summation formulas
  of Johnson, Kotz and Balakrishnan (1995),
  (31.15) on page 516 and (31.15') on p.519, the latter being typo-corrected
  for a missing  factor \eqn{1 / j!}.

  \code{dntJKBf()} is \code{\link{Vectorize}(dntJKBf1,
    c("x","df","ncp"))}, i.e., works vectorized in all three main
  arguments \code{x}, \code{df} and \code{ncp}.

  The functions \code{.dntJKBch1()} and  \code{.dntJKBch()} are only there
  for didactical reasons allowing to check that indeed formula (31.15')
  in the reference is missing a \eqn{j!} factor in the denominator.

  The \code{dntJKBf*()} functions are written to also work with
  arbitrary precise numbers of \code{\link{class}}
  \code{"\link[Rmpfr:mpfr-class]{mpfr}"} (from package \CRANpkg{Rmpfr})
  as arguments.
  %% FIXME: check that
}
\usage{
dntJKBf1(x, df, ncp, log = FALSE, M = 1000)
dntJKBf (x, df, ncp, log = FALSE, M = 1000)

## The "checking" versions, only for proving correctness of formula:
.dntJKBch1(x, df, ncp, log = FALSE, M = 1000, check=FALSE, tol.check = 1e-7)
.dntJKBch (x, df, ncp, log = FALSE, M = 1000, check=FALSE, tol.check = 1e-7)
}
\arguments{
  \item{x, df, ncp}{see \R's \code{\link{dt}()}; note that each can be
    of class \code{"\link[Rmpfr:mpfr-class]{mpfr}"}.}
  \item{log}{as in \code{\link{dt}()}, a logical indicating if
    \eqn{\log(f(x,*))} should be returned instead of \eqn{f(x,*)}.}
  \item{M}{the number of terms to be used, a positive integer.}
  \item{check}{logical indicating if checks of the formula equialities
    should be done.}
  \item{tol.check}{tolerance to be used for \code{\link{all.equal}()}
    when \code{check} is true.}
}
\details{
  How to choose \code{M} optimally has not been investigated yet.

  Note that relatedly,

  \R's source code \file{R/src/nmath/dnt.c} has claimed from 2003 till 2014
  but \bold{wrongly} that the noncentral t density \eqn{f(x, *)} is \preformatted{
    f(x, df, ncp) =
 	df^(df/2) * exp(-.5*ncp^2) /
 	(sqrt(pi)*gamma(df/2)*(df+x^2)^((df+1)/2)) *
 	   sum_{k=0}^Inf  gamma((df + k + df)/2)*ncp^k / prod(1:k)*(2*x^2/(df+x^2))^(k/2) .}

  These functions (and this help page prove that it was wrong.
}
\value{
  a number for \code{dntJKBf1()} and \code{.dntJKBch1()}.

  a numeric vector of the same length as the maximum of the lengths of
  \code{x, df, ncp} for  \code{dntJKBf()} and \code{.dntJKBch()}.
}
\references{
  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995)
  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley.
  \cr Chapter 31, Section \emph{5  Distribution Function}, p.514 ff
}
\author{Martin Maechler}
\seealso{
  \code{\link{dt}}.
}
\examples{
tt <- seq(0, 10, len = 21)
ncp <- seq(0, 6, len = 31)
dt3R   <- outer(tt, ncp, dt,     df = 3)
dt3JKB <- outer(tt, ncp, dntJKBf, df = 3)
all.equal(dt3R, dt3JKB)

% This is already in ../tests/t-nonc-tst.R --- TODO: drop here or there ?!
x <- seq(-1,12, by=1/16)
fx <- dt(x, df=3, ncp=5)
re1 <- 1 - .dntJKBch(x, df=3, ncp=5) / fx # with warnings
re2 <- 1 -  dntJKBf (x, df=3, ncp=5) / fx
stopifnot(all.equal(re1[!is.na(re1)], re2[!is.na(re1)], tol=1e-6))
matplot(x, log10(abs(cbind(re1, re2))), type = "o", cex = 1/4)
}
\keyword{distribution}
\keyword{math}
