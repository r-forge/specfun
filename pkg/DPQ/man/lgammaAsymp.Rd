\name{lgammaAsymp}
\alias{lgammaAsymp}
\title{Asymptotic Log Gamma Function}
\description{
   Compute an n-th order asymptotic approximation to log Gamma function,
   using Bernoulli numbers \code{\link{Bern}(k)} for \code{k} in
   \eqn{1, \ldots, 2n}.
}
\usage{
lgammaAsymp(x, n)
}
\arguments{
  \item{x}{numeric vector}
  \item{n}{integer specifying the approximation order.}
}
%% \details{

%% }
\value{
  numeric vector with the same attributes (\code{\link{length}()} etc) as
  \code{x}, containing approximate \code{\link{lgamma}(x)} values.
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
 \code{\link{lgamma}}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, n)
{
    s <- (x - 1/2) * log(x) - x + log(2 * pi)/2
    if (n >= 1) {
        Ix2 <- 1/(x * x)
        k <- 1:n
        Bern(2 * n)
        Bf <- rev(.bernoulliEnv$.Bern[k]/(2 * k * (2 * k - 1)))
        bsum <- Bf[1]
        for (i in k[-1]) bsum <- Bf[i] + bsum * Ix2
        s + bsum/x
    }
    else s
  }
}
\keyword{distribution}
\keyword{math}
