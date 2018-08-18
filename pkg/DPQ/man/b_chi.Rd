\name{b_chi}
\alias{b_chi}
\alias{b_chiAsymp}
\title{Compute \eqn{E[χ_ν] / \sqrt{ν}} useful for Non-central t-Distribution}
\description{
  \deqn{b_\chi(\nu) := E[\chi(\nu)] / \sqrt{\nu},}
  where \eqn{\chi(\nu)} denotes a chi-distributed random variable, i.e.,
  the square of a chi-squared variable.

  This in a relatively important auxiliary function when computing with
  non-central t distribution functions, specifically see
  Johnson, Kotz, Bala... , p.520, after (31.26a).
}
\usage{
b_chi     (nu, one.minus = FALSE, c1 = 341, c2 = 1000)
b_chiAsymp(nu, one.minus = FALSE)
}

\arguments{
  \item{nu}{non-negative numeric vector of degrees of freedom.}
  \item{one.minus}{logical indicating if \eqn{1 - b()} should be
    returned instead of \eqn{b()}.}
  \item{c1, c2}{boundaries for different approximation intervals used:
    \cr for \code{0  < nu <= c1}, internal \code{b1()} is used,
    \cr for \code{c1 < nu <= c2}, internal \code{b2()} is used, and
    \cr for \code{c2 < nu}, the \code{b_chiAsymp()} function is used,
    (and you can use that explicitly, also for smaller \code{nu}).
  }
}
\value{
  a numeric vector of the same length as \code{nu}.
}
\references{
  Johnson, Kotz,  Balakrishnan (1995)
  \emph{Continuous Univariate Distributions},
  Vol 2, 2nd Edition; Wiley.
  Formula on page 520, after (31.26a)
}
\author{Martin Maechler}
\seealso{
  The t-distribution (base \R) page \code{\link{pt}}.
}
\examples{
r <- curve(b_chi, 1e-10, 1e5, log="x")
with(r, lines(x, b_chi(x, one.minus=TRUE), col = 2))

## Zoom in to c1-region
rc1 <- curve(b_chi, 340.5, 341.5, n=1001)

## Zoom in to c2-region
rc2 <- curve(b_chi, 999.5, 1001.5, n=1001)

grCol <- adjustcolor("forest green", 1/2)
curve(b_chi,                    1/2, 1e11, log="x")
curve(b_chiAsymp, add = TRUE, col = grCol, lwd = 3)
## 1-b(nu) ~= 1/(4 nu) a power function <==> linear in log-log scale:
curve(b_chi(x, one.minus=TRUE), 1/2, 1e11, log="xy")
curve(b_chiAsymp(x, one.minus=TRUE), add = TRUE, col = grCol, lwd = 3)

%% FIXME: and see an example in package 'DPQmpfr'
}
\keyword{math}
\keyword{distribution}
