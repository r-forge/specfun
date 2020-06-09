\name{phyperIbeta}
\alias{phyperIbeta}
\title{Pearson's incomplete Beta Approximation to the Hyperbolic Distribution}
\description{
  Pearson's incomplete Beta function approximation to the cumulative
  hyperbolic distribution function \code{\link{phyper}(.)}.

  Note that in \R, \code{\link{pbeta}()} provides a version of the
  incomplete Beta function.
}
\usage{
phyperIbeta(q, m, n, k)
}
\arguments{
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in
    \eqn{0,1,\dots, m+n}.}
}

\value{
  a numeric vector \dQuote{like} \code{q} with values approximately equal
  to \code{\link{phyper}(q,m,n,k)}.
}

\references{
  Johnson, Kotz & Kemp (1992):  (6.90), p.260 -- Bol'shev (1964)
}
\author{Martin Maechler}
\seealso{
 \code{\link{phyper}}.
}
\examples{
## The function is currently defined as
function (q, m, n, k)
{
    Np <- m
    N <- n + m
    n <- k
    x <- q
    p <- Np/N
    np <- n * p
    xi <- (n + Np - 1 - 2 * np)/(N - 2)
    d.c <- (N - n) * (1 - p) + np - 1
    cc <- n * (n - 1) * p * (Np - 1)/((N - 1) * d.c)
    lam <- (N - 2)^2 * np * (N - n) * (1 - p)/((N - 1) * d.c *
        (n + Np - 1 - 2 * np))
    pbeta(1 - xi, lam - x + cc, x - cc + 1)
  }
}
\keyword{ distribution }
