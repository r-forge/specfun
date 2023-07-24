\name{phyperBin}
\title{HyperGeometric Distribution via Approximate Binomial Distribution}
\alias{phyperBin.1}
\alias{phyperBin.2}
\alias{phyperBin.3}
\alias{phyperBin.4}
\description{%% ../R/hyper-dist.R
  Compute hypergeometric cumulative probabilities via (good) binomial
  distribution approximations.
  The arguments of these functions are \emph{exactly} those of \R's own
  \code{\link{phyper}()}.
}% ~/R/D/r-devel/R/src/library/stats/man/Hypergeometric.Rd
\usage{
phyperBin.1(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBin.2(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBin.3(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
phyperBin.4(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in \eqn{0,1,\dots, m+n}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}

\details{
  TODO
}
\value{
 a \code{\link{numeric}} vector, with the length the maximum of the
 lengths of \code{q, m, n, k}.
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
  \code{\link{phyper}},
  \code{\link{pbinom}}
}
\examples{
## The 1st function is
function (q, m, n, k, lower.tail = TRUE, log.p = FALSE)
  pbinom(q, size = k, prob = m/(m + n), lower.tail = lower.tail,
         log.p = log.p)
}
\keyword{distribution}
