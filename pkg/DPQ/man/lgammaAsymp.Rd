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
\author{Martin Maechler}
\seealso{
 \code{\link{lgamma}}; the \eqn{n}-th Bernoulli number
 \code{\link{Bern}(n)}, and also \emph{exact} fractions Bernoulli numbers
 \code{\link[gmp]{BernoulliQ}()} from package \CRANpkg{gmp}.
}
\examples{
## The function is currently
lgammaAsymp

}
\keyword{distribution}
\keyword{math}
