\name{qnormR}
\title{Pure R version of \R's \code{qnorm()} with Diagnostics and Tuning Parameters}
\alias{qnormR1}
\alias{qnormR}
\description{
  Compute's \R level implementations of \R's \code{\link{qnorm}()} as
  implemented in C code (in \R's \file{mathlib}).
}
\usage{
qnormR1(p, mu = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, trace = 0)
qnormR (p, mu = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, trace = 0)
}
\arguments{
  \item{p}{probability \eqn{p}, \eqn{1-p}, or \eqn{\log(p)},
    \eqn{\log(1-p)}, depending on \code{lower.tail} and \code{log.p}.}
  \item{mu}{mean of the normal distribution.}
  \item{sd}{standard deviation of the normal distribution.}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qnorm}()}.}
  \item{trace}{logical or integer; if positive or \code{TRUE}, diagnostic
    output is printed to the console during the computations.}
}
\details{
  For \code{qnormR1(p, ..)}, \code{p} must be of length one, whereas
  \code{qnormR(p, m, s, ..)} works vectorized in \code{p}, \code{mu}, and
  \code{sd}.  In the \pkg{DPQ} package source, it simply the result of
  \code{\link{Vectorize}(qnormR1, ...)}.
}
\value{
  a numeric vector like the input \code{q}.
}
%% \references{
%% }
\author{Martin Maechler}
%% \note{
%% }
\seealso{
  \code{\link{qnorm}}
}
\examples{
qR <- curve(qnormR, n = 2^11)
abline(h=0, v=0:1, lty=3, col=adjustcolor(1, 1/2))
with(qR, all.equal(y, qnorm(x), tol=0)) # currently shows TRUE
with(qR, all.equal(pnorm(y), x, tol=0)) # currently: mean rel. diff.: 2e-16
stopifnot(with(qR, all.equal(pnorm(y), x, tol = 1e-14)))
}
\keyword{distribution}
