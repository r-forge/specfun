\name{qtR}
\alias{qtR}
\alias{qtR1}
\title{Pure \R Implementation of \R's C-level t-Distribution Quantiles \code{qt()}}
\description{
  A pure \R implementation of \R's Mathlib own C-level \code{\link{qt}()} function.
  \cr % in  ../R/t-dist.R
  \code{qtR()} is simply defined as \preformatted{qtR <- Vectorize(qtR1, c("p","df"))}
  where in \code{qtR1(p, df, *)} both \code{p} and \code{df} must be of length one.
}
\usage{
qtR1(p, df, lower.tail = TRUE, log.p = FALSE, .....)
qtR (p, df, lower.tail = TRUE, log.p = FALSE, .....)
}
\arguments{
  \item{p, df}{vectors of probabilities and degrees of freedom, see \code{\link{qt}}.}
  \item{lower.tail, log.p}{logical; see \code{\link{qt}}.}
  %% \item{tol}{non-negative convergence tolerance passed to \code{\link{uniroot}()}.}
  %% \item{verbose}{logical indicating if \emph{every} call of the objective
  %%   function should produce a line of console output.}
  %% \item{\dots}{optional further arguments passed to \code{\link{uniroot}()}.}
}
\value{
  numeric vector of t quantiles, properly recycled in \code{(p, df)}.
}
\author{Martin Maechler}
\seealso{
 \code{\link{qtU}} and \R's \code{\link{qt}}.
}
\examples{
## Inspired from Bugzilla PR#16380
pxy <- curve(pt(-x, df = 1.09, log.p = TRUE), 4e152, 1e156, log="x", n=501)
qxy <- curve(-qt(x, df = 1.09, log.p = TRUE), -392, -385, n=501, log="y", col=4, lwd=2)
lines(x ~ y, data=pxy, col = adjustcolor(2, 1/2), lwd=5, lty=3)
## now our "pure R" version:
qRy <- -qtR(qxy$x, df = 1.09, log.p = TRUE)
all.equal(qRy, qxy$y)

try( # fails
qtR1(.1, .1, verbose=TRUE)
)
}
