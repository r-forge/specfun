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
}
\keyword{distribution}
\keyword{math}
