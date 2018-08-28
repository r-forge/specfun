\name{lgamma1p}
\title{Accurate \code{log(gamma(a+1))}}
\alias{lgamma1p}
\alias{lgamma1p.}
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
lgamma1p (a, tol_logcf = 1e-14)     % ../R/beta-fns.R
lgamma1p.(a, cutoff.a = 1e-6, k = 3)% ../R/qchisqAppr.R
lgamma1p_series(x, k)% ../R/qchisqAppr.R
}
\arguments{
  \item{a, x}{a numeric vector.}
  \item{tol_logcf}{for \code{lgamma1p()}: a non-negative number ...}% FIXME
  \item{cutoff.a}{for \code{lgamma1p.()}: a positive number indicating
    the cutoff to switch from  ...}% FIXME
  \item{k}{an integer, the number of terms in the series expansion used internally.}
}
\details{
  \code{lgamma1p()} is an \R translation of the function (in Fortran) in
  Didonato and Morris (1992) which uses a 40-degree polynomial approximation.

  \code{lgamma1p_series(x, k)} is Taylor series approximation of order \code{k},
  (derived via Maple), which is \eqn{-\gamma x + \pi^2 x^2/ 12 +
  O(x^3)}{gamma*x + pi^2 * x^2/ 12 + O(x^3)}, where \eqn{\gamma}{gamma}
  is Euler's constant 0.5772156649....
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
\author{Martin Maechler, for \code{lgamma1p()}, translated from the
  C/Fortran code of algorithm 708 by Didonato and Morris.
}
%% \note{
%% }
\seealso{
  \code{\link{log1p}}, \code{\link{pbeta}}.
}
\examples{
curve(-log(x*gamma(x)), 1e-30, .8, log="xy", col="gray50", lwd = 3,
      axes = FALSE, ylim = c(1e-30,1))
sfsmisc::eaxis(1); sfsmisc::eaxis(2)
at <- 10^(1-4*(0:8))
abline(h = at, v = at, col = "lightgray", lty = "dotted")
curve(-lgamma( 1+x), add=TRUE, col="red2", lwd=1/2)# underflows even earlier
curve(-lgamma1p (x), add=TRUE, col="blue")
curve(-lgamma1p.(x), add=TRUE, col=adjustcolor("forest green",1/4),
      lwd = 5, lty = 2)
for(k in 1:6)
 curve(-lgamma1p_series(x, k=k), add=TRUE, col=paste0("gray",20+k*10), lty = 3)
}
\keyword{distribution}
\keyword{math}

