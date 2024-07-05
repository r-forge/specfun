\name{Ixpq}
\alias{Ixpq}
\title{Normalized Incomplete Beta Function "Like" \code{pbeta()}}
\description{
  Computes the normalized incomplete beta function, in pure \R code,
  derived from Nico Temme's Maple code for computing Table 1 in  Gil et al (2023).

  It uses a continued fraction, similarly to \code{bfrac()} in the TOMS 708
  algorithm underlying \R's  \code{\link{pbeta}()}.
}
\usage{
Ixpq(x, l_x, p, q, tol = 3e-16, it.max = 100L, plotIt = FALSE)
}
\arguments{
  \item{x}{numeric}%% FIXME should work for mpfr-numbers too !!
  \item{l_x}{\code{1 - x}; may be specified with higher precision (e.g.,
    when \eqn{x \approx 1}, \eqn{1-x} suffers from cancellation).}
  \item{p, q}{the two shape parameters of the beta distribution.}
  \item{tol}{positive number, the convergence tolerance for the continued fraction computation.}
  \item{it.max}{maximal number of continued fraction steps.}
  \item{plotIt}{a \code{\link{logical}}, if true, plots show the relative
    approximation errors in each step.}
}
%% \details{
%% }
\value{
  a vector like \code{x} or \code{l_x} with corresponding
  \code{\link{pbeta}(x, *)} values.
}
\references{
  Gil et al. (2023) %% FIXME add --> ../Misc/pnt-Gil_etal-2023/
}
\author{Martin Maechler; based on original Maple code by Nico Temme.}
%% \note{
%% %%  ~~further notes~~
%% }
\seealso{
 \code{\link{pbeta}}, \code{\link{pbetaRv1}()}, ..
}
\examples{
x <- seq(0, 1, by=1/16)
r <- Ixpq(x, 1-x, p = 4, q = 7, plotIt = TRUE)
cbind(x, r)
## and "test" ___FIXME__
}
\keyword{math}
\keyword{distribution}
