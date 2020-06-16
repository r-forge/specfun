\name{betaD94}
\title{Ding(1994) (non-central) Beta Distribution Functions}
\alias{dbetaD94}
\alias{pbetaD94}
\alias{qbetaD94}
\description{
  The three cases use Ding(1994)'s algorithm A, B, and C, respectively,
  each of which implements a recursion formula using only simple
  arithmetic and \code{\link{log}} and \code{\link{exp}}.

  These are particularly useful also for using with high precision
  \code{"mpfr"} numbers from the \CRANpkg{Rmpfr} CRAN package.

>> ../R/pbeta-Ding94.R -- including examples (for  examples /
   -------------------     ../tests/ )
}
\usage{
dbetaD94(x, shape1, shape2, ncp = 0, log = FALSE,
         eps = 1e-10, itrmax = 100000L, verbose = TRUE)
pbetaD94(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
         eps = 1e-10, itrmax = 100000L, verbose = TRUE)
qbetaD94(p, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
         delta = 1e-6,
         eps = delta^2,
         itrmax = 100000L,
         iterN = 1000L,
         verbose = TRUE)
}

\arguments{
  \item{x,q}{numeric vector of values in \eqn{[0,1]} as beta variates.}
  \item{shape1, shape2}{the two shape parameters of the beta
    distribution, must be positive.}
  \item{ncp}{the noncentrality parameter; by default zero for the
    (\emph{central}) beta distribution; if positive, we have a noncentral
    beta distribution.}
  \item{p}{numeric vector of probabilities, \code{\link{log}()}ged in
    case \code{log.p} is true.}
  \item{log, log.p}{logical indicating if the density or probability
    values should be \code{\link{log}()}ged.}
  \item{lower.tail}{logical indicating if the lower or upper tail
    probability should be computed, or for \code{qbeta*()} are provided.}
  \item{eps}{a non-negative number specifying the desired accuracy for
    computing F() and f().}
  \item{itrmax}{the maximal number of steps for computing F() and f().}
  \item{delta}{[For \code{qbeta*()}:] non-negative number indicating the
    desired accuracy for computing \eqn{x_p} (the root of \eqn{pbeta*() == p}),
    i.e., the convergence tolerance for the Newton iterations.}
  \item{iterN}{[For \code{qbeta*()}:] The maximal number of \bold{N}ewton iterations.}
  \item{verbose}{logical (or integer) indicating the amount of
    diagnostic output during computation.}
}
%% \details{
%%
%% }
\value{
  In all three cases, a numeric vector with the same attributes as
  \code{x},
  containing (an approximation) to the correponding beta distribution function.
}
\references{
  Cherng G. Ding (1994) % ~/save/papers/Numerics/Ding_1994_noncentr-beta.pdf
  On the computaiton of the noncentral beta distribution.
  Computational Statistics & Data Analysis \bold{18}, 449--455.
}
\author{Martin Maechler}
\seealso{
  \code{\link{pbeta}}.
}
\examples{
 " TODO FIXME "
}
\keyword{distribution}
\keyword{math}
