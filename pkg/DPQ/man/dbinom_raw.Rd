\name{dbinom_raw}
\title{R's C (Mathlib) dbinom_raw() Binomial Probability Function in pure R}
\alias{dbinom_raw}% >> ../R/dbinom-nbinom.R
\description{
  A pure \R implementation R's C API (\sQuote{Mathlib} specifically)
  \code{dbinom_raw()} function which computes binomial probabilities
  \emph{and} is continuous in \code{x}, i.e., also \dQuote{works} for
  non-integer \code{x}.
}
\usage{
dbinom_raw (x, n, p, q, log = FALSE)
}
\arguments{
  \item{x}{vector with values typically in \code{0:n},  but here allowed to
    non-integer values.}
  \item{n}{called \code{size} in \R's \code{\link{dbinom}()}.}
  \item{p}{called \code{prob} in \R's \code{\link{dbinom}()}, the success
    probability, hence in \eqn{[0, 1]}.}
  \item{q}{mathemtically the same as \eqn{1 - p}, but may be (much) more
    accurate, notably when small.}
  \item{log}{logical indicating if the \code{\link{log}()} of the resulting
    probability should be returned; useful notably in case the probability
    itself would underflow to zero.}
}
\value{
  numeric vector of the same length as \code{x} (which may have to be
  thought of recycled along \code{n}, \code{p} and/or \code{q}.
}
\author{R Core and Martin Maechler}
\seealso{
  Note that our CRAN package \CRANpkg{Rmpfr} provides
  \code{\link[Rmpfr]{dbinom}}, an mpfr-accurate function to be used
  used instead of \R's or this pure \R version relying \code{\link{bd0}()} and
  \code{\link{stirlerr}()} where the latter currently only provides
  accurate double precision accuracy.
}
\examples{%% for now have regr.tests only in ../tests/dnbinom-tst.R
}
\keyword{distribution}
