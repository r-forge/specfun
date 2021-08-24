\name{log1pmx}
\alias{log1pmx}
\alias{log1pmxC}
\title{Accurate  \code{log(1+x) - x}  Computation}
\description{
  Compute \deqn{\log(1+x) - x}{log(1+x) - x}
  accurately also for small \eqn{x}, i.e., \eqn{|x| \ll 1}{|x| << 1}.

  Since April 2021, the pure \R code version \code{log1pmx()} also works
  for "mpfr" numbers (from package \CRANpkg{Rmpfr}).
}
\usage{
log1pmx (x, tol_logcf = 1e-14, eps2 = 0.01, minL1 = -0.79149064,
         trace.lcf = FALSE,
         logCF = if(is.numeric(x)) logcf else logcfR.)
log1pmxC(x)  # TODO in future: arguments (minL1, eps2, tol_logcf),
             # possibly with *different* defaults (!)
}
\arguments{
  \item{x}{numeric (or \code{"\link[Rmpfr]{mpfr}"} number) vector with values \eqn{x > -1}.}
  \item{tol_logcf}{a non-negative number indicating the tolerance
    (maximal relative error) for the auxiliary \code{\link{logcf}()} function.}
  \item{eps2}{positive cutoff where the algorithm switches from a few
    terms, to using \code{\link{logcf}()} explicitly.  Note that for
    more accurate mpfr-numbers the default \code{eps = .01} is too large,
    even more though when the tolerance is lowered (from \code{1e-14}).}
  \item{minL1}{negative cutoff, called \code{minLog1Value} in Morten
    Welinder's C code for \code{log1pmx()} in \file{R/src/nmath/pgamma.c},
    hard coded there to -0.79149064 which seems not optimal for
    computation of \code{log1pmx()}, at least in some cases, and hence
    \bold{the default may be changed in the future}.}
  \item{trace.lcf}{\code{\link{logical}} used in
    \code{\link{logcf}(.., trace=trace.lcf)}.}
  \item{logCF}{the \code{\link{function}} to be used as
    \code{\link{logcf}()}.  The default chooses the pure \R \code{logcfR()}
    when \code{x} is not numeric, and chooses the C-based \code{logcf()}
    when \code{is.numeric(x)} is true.}
}
\details{% ../R/beta-fns.R
  In order to provide full accuracy,
  the computations happens differently in three regions for \eqn{x},
  \deqn{m_l = \code{minL1} = -0.79149064} is the first cutpoint,
  \describe{
    \item{\eqn{x < m_l} or \eqn{x > 1}:}{use \code{log1pmx(x) := \link{log1p}(x) - x},}
    \item{\eqn{|x| < \epsilon_2}{|x| < eps2}:}{use
      \eqn{t((((2/9 * y + 2/7)y + 2/5)y + 2/3)y - x)},}
    \item{\eqn{x \in [ml,1]}, and \eqn{|x| \ge \epsilon_2}{|x| >= eps2}:}{use
      \eqn{t(2y logcf(y, 3, 2) - x)},}
  }
  where \eqn{t := \frac{x}{2 + x}}{t := x/(2+x)}, and \eqn{y := t^2}.

  Note that the formulas based on \eqn{t} are based on the (fast
  converging) formula
  \deqn{\log(1+x) = 2\left(r + \frac{r^3}{3}+ \frac{r^5}{5} + \ldots\right),}{%
    log(1+x) = 2(r + r^3/3  + r^5/5 + ...),}
  where \eqn{r := x/(x+2)}, see the reference.

  \code{log1pmxC()} is an interface to \R C API (\file{Rmathlib}) function.
}
\value{a numeric vector (with the same attributes as \code{x}).
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.\cr
  Formula (4.1.29), p.68.

  Martin \enc{Mächler}{Maechler} (2021).
  %log1pmx(), bd0(), stirlerr() -- Computing Poisson, Binomial, Gamma Probabilities in \R}
  \emph{log1pmx, \dots Computing \dots Probabilities in \R}.
  (\pkg{DPQ} package vignette)% ../vignettes/log1pmx-etc.Rnw
  % FIXME -- not yet:
  % \url{https://CRAN.R-project.org/package=DPQ/vignettes/log1pmx-etc.pdf}.
}
\author{A translation of Morten Welinder's C code of Jan 2005, see R's bug
  issue \PR{7307}, parametrized and tuned by Martin Maechler.
}
\seealso{
  \code{\link{logcf}}, the auxiliary function,
  \code{\link{lgamma1p}} which calls \code{log1pmx}, \code{\link{log1p}}
}
\examples{
(doExtras <- DPQ:::doExtras()) # TRUE e.g. if interactive()
n1 <- if(doExtras) 1001 else 201
curve(log1pmx, -.9999, 7, n=n1); abline(h=0, v=-1:0, lty=3)
curve(log1pmx, -.1,  .1,  n=n1); abline(h=0, v=0, lty=3)
curve(log1pmx, -.01, .01, n=n1) -> l1xz2; abline(h=0, v=0, lty=3)
## C and R versions correspond closely:
with(l1xz2, stopifnot(all.equal(y, log1pmxC(x), tol = 1e-15)))

e <- if(doExtras) 2^-12 else 2^-8; by.p <- 1/(if(doExtras) 256 else 64)
xd <- c(seq(-1+e, 0+100*e, by=e), seq(by.p, 5, by=by.p)) # length 676 or 5476 if do.X.
plot(xd, log1pmx(xd), type="l", col=2, main = "log1pmx(x)")
abline(h=0, v=-1:0, lty=3)

## much more graphics etc in ../tests/dnbinom-tst.R  (and the vignette, see above)
}
\keyword{math}
