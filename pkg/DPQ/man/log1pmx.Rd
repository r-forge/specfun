\name{log1pmx}
\alias{log1pmx}
\alias{log1pmxC}
\title{Accurate \code{log(1+x) - x}}
\description{
  Compute \deqn{\log(1+x) - x}{log(1+x) - x}
  accurately also for small \eqn{x}, i.e., \eqn{|x| \ll 1}{|x| << 1}.

  Since April 2021, the pure \R code version \code{log1pmx()} also works
  for "mpfr" numbers (from package \CRANpkg{Rmpfr}).
}
\usage{
log1pmx (x, tol_logcf = 1e-14, eps2 = 0.01, minL1 = -0.79149064,
         trace.lcf = FALSE,
         logCF = if(is.numeric(x)) logcf else logcfR)
log1pmxC(x)
}
\arguments{
  \item{x}{numeric vector with values \eqn{x > -1}.}
  \item{tol_logcf}{a non-negative number indicating the tolerance
    (maximal relative error) for the auxiliary \code{\link{logcf}()} function.}
  \item{eps2}{positive cutoff where the algorithm switches from a few
    terms, to using \code{\link{logcf}()} explicitly.}
  \item{minL1}{negative cutoff, called \code{minLog1Value} in Morten
    Welinder's C code for \code{log1pmx()} in \file{R/src/nmath/pgamma.c}.
    It seems to Martin M that it is not quite optimal, at least in some cases.}
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
}
\author{A translation of Morten Welinder's C code of Jan 2005, see R's bug
  issue \PR{7307}.
}
\seealso{
  \code{\link{logcf}}, the auxiliary function,
  \code{\link{lgamma1p}} which calls \code{log1pmx}, \code{\link{log1p}}
}
\examples{
(doExtras <- DPQ:::doExtras()) # TRUE e.g. if interactive()
n1 <- if(doExtras) 1001 else 201
l1x <- curve(log1pmx, -.9999, 7, n=n1)
abline(h=0, v=-1:0, lty=3)
l1xz  <- curve(log1pmx, -.1,  .1,  n=n1); abline(h=0, v=0, lty=3)
l1xz2 <- curve(log1pmx, -.01, .01, n=n1); abline(h=0, v=0, lty=3)
l1xz3 <- curve(-log1pmx(x), -.002, .002, n=n1, log="y", yaxt="n")
sfsmisc::eaxis(2); abline(v=0, lty=3)
with(l1xz3, stopifnot(all.equal(y, -log1pmxC(x))))

e <- if(doExtras) 2^-12 else 2^-8; by.p <- 1/(if(doExtras) 256 else 64)
xd <- c(seq(-1+e, 0+100*e, by=e), seq(by.p, 5, by=by.p))
plot(xd, log1pmx(xd), type="l", col=2, main = "log1pmx(x)")
abline(h=0, v=-1:0, lty=3)

if(requireNamespace("Rmpfr") && packageVersion("sfsmisc") >= "1.1-10") withAutoprint({
  xM <- Rmpfr::mpfr(xd, 512)
  ## for MPFR numbers, really need more than tol_logcf = eps = 1e-14 (default)
  system.time(lg1pM <- log1pmx(xM, tol_logcf = 1e-25, eps2 = 1e-4))
              # 4.5 sec if(doExtras) 0.43s otherwise -- *NOT* really needed!
  ## MM: But really, this should be more accurate anyway (?!):
  lg1pM. <- log1p(xM) - xM
  xM2k <- Rmpfr::mpfr(xd, 2048)
  lg1pM2k <- log1p(xM2k) - xM2k # even more accurate
  relErrV <- sfsmisc::relErrV
  asNumeric <- Rmpfr::asNumeric
  reE00 <- asNumeric( relErrV(lg1pM2k, lg1pM.) )
  print(signif(range(reE00)), 2) # [-1.5e-151, 4.8e-151] -- 512 bits is perfect
  ## the error of the log1pmx() "algorithm" even for "perfect" mpfr-accuracy:
  rE.log1pm <- asNumeric(relErrV(lg1pM2k, lg1pM))
  print(signif(range(  rE.log1pm ), 2)) # -1.2e-27 3.7e-49

  re <- asNumeric(relErrV(lg1pM., log1pmx(xd)))
  plot(xd, re, type="b", cex=1/2)
  abline(h = (-2:2)*2^-52, lty=2, col=adjustcolor("gray20", 1/2))
  ## only negative x:
  iN <- xd < 0
  iN <- -0.84 < xd & xd < -0.4
  plot(xd[iN], re[iN], type="b", cex=1/2)
  abline(h = (-2:2)*2^-52, lty=2, col=adjustcolor("gray20", 1/2))

  plot(xd[iN], abs(re[iN]), type="b", cex=1/2, log="y",
       main = "| relErr( log1pmx(x) ) |  {via 'Rmpfr'}",
       ylim = c(4e-17, max(abs(asNumeric(re)[iN]))))
  abline(h = c(1:2,4)*2^-52, lty=2, col=adjustcolor("gray20", 1/2))
  mL1 <- eval(formals(log1pmx)$minL1)
  abline(v = mL1, lwd=3, col=adjustcolor(2, 1/2))
  axis(3, at=mL1, "minL1", col.axis=2, col=2)

  re2 <- asNumeric(sfsmisc::relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = -0.7)))
  lines(xd[iN], abs(re2), col=adjustcolor(4, 1/2), lwd=2)
  abline(v = -0.7, lwd=3, col=adjustcolor(4, 2/3), lty=3)
  axis(3, line=-1, at=-0.7, "mL1 = -0.7", col.axis=4, col=4) ## it seems that would be better

  re3 <- asNumeric(sfsmisc::relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = -0.67)))
  lines(xd[iN], abs(re3), col=adjustcolor(6, 1/3), lwd=2)
  abline(v = -0.67, lwd=3, col=adjustcolor(6, 2/3), lty=3)
  axis(3, line=-2, at=-0.67, "mL1 = -0.67", col.axis=6, col=6)
  lines(lowess(xd[iN], abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
   # MM: I'm confused -- why does -0.7 not show problems to the right, but
   #                     but      -0.67 *does* show them .. ?


  re4 <- asNumeric(sfsmisc::relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = -0.64)))
  lines(xd[iN], abs(re4), col=adjustcolor(7, 1/3), lwd=2)
  abline(v = -0.64, lwd=3, col=adjustcolor(7, 2/3), lty=3)
  axis(3, line=-2, at=-0.64, "mL1 = -0.64", col.axis=7, col=7)
  lines(lowess(xd[iN], abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re4),    f=1/50), col=adjustcolor(7, 1/2), lwd=6)
                             # ...

  re5 <- asNumeric(sfsmisc::relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = -0.6)))
  lines(xd[iN], abs(re5), col=adjustcolor(8, 1/3), lwd=2)
  abline(v = -0.6, lwd=3, col=adjustcolor(8, 2/3), lty=3)
  axis(3, line=-1, at=-0.6, "mL1 = -0.6", col.axis=8, col=8)
  lines(lowess(xd[iN], abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re4),    f=1/50), col=adjustcolor(7, 1/2), lwd=6)
  lines(lowess(xd[iN], abs(re5),    f=1/50), col=adjustcolor(8, 1/2), lwd=6)

   # -0.6  now is clearly too large, -0.7 was better

})# if "Rmpfr"

}
\keyword{math}
