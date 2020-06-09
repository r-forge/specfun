\name{phyperR}
\alias{phyperR}
\title{\R-only version of \R's original phyper() algorithm}
\description{
  An \R version of the first \code{phyper()} algorithm in \R, which was
  used up to svn rev \code{30227} on 2004-07-09.
}
\usage{% ../R/hyper-dist.R (at end)
phyperR(q, m, n, k)
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
%% \details{

%% }
\note{
  The original argument list in \code{C} was \code{(x, NR, NB, n)} where
  there were \emph{red} and \emph{black} balls in the urn.

  Note that we have \emph{vectorized} a translation to \R of the original C
  code.
}
\value{
  a numeric vector similar to \code{phyper(q, m, n, k)}.
}
%% \references{
%% }
\author{Martin Maechler}
\seealso{
 \code{\link{phyper}}.
}
\examples{
m <- 9:12; n <- 7:10; k <- 10
x <- 0:(k+1) # length 12
for(. in 0:11) stopifnot(
  all.equal(phyper (., m=10, n=8, k=10),
            phyperR(., m=10, n=8, k=10), tol=1e-12))
## BUT the vectorization fails badly : FIXME! %% <<<<<< FIXME <<<<<<<<<<<<<<<
(A <- cbind(x, m, n, k)) # shows the recycling
cbind(A, ph = phyper (x, m, n, k),
      phR   = phyperR(x, m, n, k))
}
\keyword{ distribution }
