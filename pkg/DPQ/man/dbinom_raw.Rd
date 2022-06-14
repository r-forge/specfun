\name{dbinom_raw}
\title{R's C Mathlib (Rmath) dbinom_raw() Binomial Probability pure R Function}
\alias{dbinom_raw}% >> ../R/dbinom-nbinom.R
\description{
  A pure \R implementation of R's C API (\sQuote{Mathlib} specifically)
  \code{dbinom_raw()} function which computes binomial probabilities
  \emph{and} is continuous in \code{x}, i.e., also \dQuote{works} for
  non-integer \code{x}.
}
\usage{% NB:  dpois_raw()	--->>> ./dgamma-utils.Rd <<<
dbinom_raw (x, n, p, q = 1-p, log = FALSE, verbose = getOption("verbose"))
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
  \item{verbose}{integer indicating the amount of verbosity of
    diagnostic output, \code{0} means no output, \code{1} more, etc.}
}
\value{
  numeric vector of the same length as \code{x} which may have to be
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
\examples{%% NB: regr.tests in >>> ../tests/dnbinom-tst.R <<<

for(n in c(3, 10, 27, 100, 500, 2000, 5000, 1e4, 1e7, 1e10)) {
 x <- if(n <= 2000) 0:n else round(seq(0, n, length.out=2000))
 p <- 3/4
 stopifnot(all.equal(dbinom_raw(x, n, p, q=1-p),
                     dbinom    (x, n, p), tol = 1e-14))
}
n <- 1024 ; x <- 0:n
plot(x, dbinom_raw(x, n, p, q=1-p) - dbinom(x, n, p), type="l", main = "|db_r(x) - db(x)|")
plot(x, dbinom_raw(x, n, p, q=1-p) / dbinom(x, n, p) - 1, type="b", log="y",
     main = "rel.err.  |db_r(x / db(x) - 1)|")
}
\keyword{distribution}
