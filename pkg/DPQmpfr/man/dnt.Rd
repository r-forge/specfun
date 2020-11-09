\name{dnt}
\title{Non-central t-Distribution Density}
\alias{dntJKBm}
%% \alias{dntWVm}
\description{
  \code{dntJKBm} is a fully \CRANpkg{Rmpfr}-ified vectorized version of
 %%\link[DPQ]{dntJKBf}()} -- logically correct, but wrong by the "link to file here" rule
  \code{\link{dntJKBf}()} from \CRANpkg{DPQ} which implements the
  summation formulas of Johnson, Kotz and Balakrishnan (1995),
  (31.15) on page 516 and (31.15') on p.519, the latter being typo-corrected
  for a missing  factor \eqn{1 / j!}.
}
\usage{
dntJKBm(x, df, ncp, log = FALSE, M = 1000)
%% dntWVm (x, df, ncp, log = FALSE)
}
\arguments{
  \item{x, df, ncp}{see \R's \code{\link{dt}()}; note that each can be
    of class \code{"\link[Rmpfr:mpfr-class]{mpfr}"}.}
  \item{log}{as in \code{\link{dt}()}, a logical indicating if
    \eqn{\log(f(x,*))} should be returned instead of \eqn{f(x,*)}.}
  \item{M}{the number of terms to be used, a positive integer.}
}
\details{
  How to choose \code{M} optimally has not been investigated yet and
  is probably also a function of the precision of the first three arguments (see
  %% \code{\link[Rmpfr]{getPrec}}. -- logically correct, but ...
  \code{\link{getPrec}} from \CRANpkg{Rmpfr}).
}
\value{
  an \code{\linkS4class{mpfr}} vector of the same length as the maximum
  of the lengths of \code{x, df, ncp}.
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
require(Rmpfr)
% FIXME: Vectorize() is not optimal, as typically df & ncp are scalar (here only df is)
## [not too large, as dntJKBm() is currently somewhat slow]
(mt  <- mpfr(tt  <- seq(0, 9, by =  1 ), 128))
(mcp <- mpfr(ncp <- seq(0, 5, by = 1/2), 128))
dt3R <- outer(tt, ncp, dt,   df = 3)
dt3M <- outer(mt, mcp, dntJKBm, df = 3, M = 128)# for speed
%% dt3V <- outer(mt, mcp, dntWVm , df = 3)
all.equal(dt3R, dt3M) # TRUE, and show difference
all.equal(dt3R, dt3M, tol=0) # 1.2e-12
%% all.equal(dt3R, dt3V, tol=0) # 1.2e-12
}
\keyword{distribution}
\keyword{math}
