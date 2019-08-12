\name{log1pmx}
\alias{log1pmx}
\title{Accurate \code{log(1+x) - x}}
\description{
  Compute \deqn{\log(1+x) - x}{log(1+x) - x}
  accurately also for small \eqn{x}, i.e., \eqn{|x| \ll 1}{|x| << 1}.
}
\usage{
log1pmx(x, tol_logcf = 1e-14)
}
\arguments{
  \item{x}{numeric vector with values \eqn{x > -1}.}
  \item{tol_logcf}{a non-negative number indicating the tolerance
    (maximal relative error) for the auxiliary \code{\link{logcf}()} function.}
}
\details{% ../R/beta-fns.R
  In order to provide full accuracy,
  the computations happens differently in three regions for \eqn{x},
  \deqn{m_l = -0.79149064} is the first cutpoint,
  \describe{
    \item{\eqn{x < ml} or \eqn{x > 1}:}{use \code{log1pmx(x) := \link{log1p}(x) - x},}
    \item{\eqn{|x| < 0.01}:}{use \eqn{t((((2/9 * y + 2/7)y + 2/5)y + 2/3)y - x)},}
    \item{\eqn{x \in [ml,1]}, and \eqn{|x| >= 0.01}:}{use \eqn{t(2y logcf(y, 3, 2) - x)},}
  }
  where \eqn{t := \frac{x}{2 + x}}{t := x/(2+x)}, and \eqn{y := t^2}.

  Note that the formulas based on \eqn{t} are based on the (fast
  converging) formula
  \deqn{\log(1+x) = 2\left(r + \frac{r^3}{3}+ \frac{r^5}{5} + \ldots\right),}{%
    log(1+x) = 2(r + r^3/3  + r^5/5 + ...),}
  where \eqn{r := x/(x+2)}, see the reference.
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
l1x <- curve(log1pmx, -.9999, 7, n=1001)
abline(h=0, v=-1:0, lty=3)
l1xz  <- curve(log1pmx, -.1, .1, n=1001); abline(h=0, v=0, lty=3)
l1xz2 <- curve(log1pmx, -.01, .01, n=1001); abline(h=0, v=0, lty=3)
l1xz3 <- curve(-log1pmx(x), -.002, .002, n=2001, log="y", yaxt="n")
sfsmisc::eaxis(2); abline(v=0, lty=3)
}
\keyword{math}
