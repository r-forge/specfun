\name{dnbinomR}
\alias{dnbinomR}% >> ../R/dbinom-nbinom.R
\alias{dnbinom.mu}
\title{Pure R Versions of R's C (Mathlib) dnbinom() Negative Binomial Probabilities}
\description{
  Compute pure \R implementations of \R's C (Mathlib)
  \code{\link{dnbinom}()} binomial probabilities, allowing to see the
  effect of the cutoff \code{eps}.
}
\usage{
dnbinomR  (x, size, prob, log = FALSE, eps = 1e-10)
dnbinom.mu(x, size, mu,   log = FALSE, eps = 1e-10)
}
\arguments{
  \item{x, size, prob, mu, log}{see \R's \code{\link[stats]{dnbinom}()}.}
  \item{eps}{non-negative number specifying the cutoff for \dQuote{small
      \code{x/size}}, in which case the 2-term approximation from
     Abramowitz and Stegun, 6.1.47 (p.257) is preferable to the
     \code{\link{dbinom}()} based evaluation.}
}
\value{
  numeric vector of the same length as \code{x} (which may have to be
  thought of recycled along \code{size} and \code{prob} or \code{mu}.
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.
}
\author{R Core and Martin Maechler}
\seealso{
  \code{\link{dbinom_raw}};
  Note that our CRAN package \CRANpkg{Rmpfr} provides
  \code{\link[Rmpfr]{dnbinom}},
  \code{\link[Rmpfr]{dbinom}} and more, where mpfr-accurate functions are
  used instead of \R's (and our pure \R version of) \code{\link{bd0}()} and
  \code{\link{stirlerr}()}.
}
\examples{%% regr.tests are in ../tests/dnbinom-tst.R
 stopifnot( dnbinomR(0, 1, 1) == 1 )
 size <- 1000 ; x <- 0:size
 dnb <- dnbinomR(x, size, prob = 5/8, log = FALSE, eps = 1e-10)
 plot(x, dnb, type="b")
 all.equal(dnb, dnbinom(x, size, prob = 5/8)) ## mean rel. diff: 0.00017...

 dnbm <- dnbinom.mu(x, size, mu = 123, eps = 1e-10)
 all.equal(dnbm, dnbinom(x, size, mu = 123)) #  Mean relative diff: 0.00069...
}
\keyword{distribution}
