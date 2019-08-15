\name{dchisqApprox}
\title{Approximations of the (Noncentral) Chisquared Density}
\alias{dnchisqR}
\alias{dchisqAsym}
\alias{dnchisqBessel}
\alias{dnoncentchisq}
\description{
  Compute the density function \eqn{f(x, *)} of the (noncentral) chisquared
  distribution.
}
\usage{%- ../R/dnchisq-fn.R
dnchisqR     (x, df, ncp, log = FALSE,
              eps = 5e-15, termSml = 1e-10, ncpLarge = 1000)
dnchisqBessel(x, df, ncp, log = FALSE)
dchisqAsym   (x, df, ncp, log = FALSE)
dnoncentchisq(x, df, ncp, kmax = floor(ncp/2 + 5 * (ncp/2)^0.5))
}
\arguments{
  \item{x}{non-negative numeric vector.}
  \item{df}{degrees of freedom (parameter), a positive number.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{log}{logical indicating if the result is desired on the log scale.}
  \item{eps}{positive convergence tolerance for the series expansion: Terms
    are added while \code{term * q > (1-q)*eps}, where \code{q} is the term's
    multiplication factor.}
  \item{termSml}{positive tolerance: in the series expansion, terms are
    added to the sum as long as they are not smaller than \code{termSml *
      sum} even when convergence according to \code{eps} had occured.  This
    was not part of the original C code, but was added later for
    safeguarding against infinite loops, from \PR{14105}, e.g., for
    \code{dchisq(2000, 2, 1000)}.}
  \item{ncpLarge}{in the case where \code{mid} underflows to \code{0}, when
    \code{log} is true, or\code{ncp >= ncpLarge}, use a central
    approximation.  In theory, an optimal choice of \code{ncpLarge} would
    not be arbitrarily set at \code{1000} (hardwired in \R's
    \code{\link{dchisq}()} here), but possibly also depend on \code{x} or
    \code{df}.}
  \item{kmax}{the number of terms in the sum for \code{dnoncentchisq()}.}
}
\details{
  \code{dnchisqR()} is a pure \R implementation of \R's own C implementation
  in the sources, \file{R/src/nmath/dnchisq.c}, additionally exposing the
  three \dQuote{tuning parameters} \code{eps}, \code{termSml}, and \code{ncpLarge}.

  \code{dnchisqBessel()} implements Fisher(1928)'s exact closed form formula
  based on the Bessel function \eqn{I_{nu}}, i.e., \R's
  \code{\link{besselI}()} function;
  specifically formula (29.4) in Johnson et al. (1995).

  \code{dchisqAsym()} is the simple asymptotic approximation from
  Abramowitz and Stegun's formula \code{26.4.27}, p. 942.

  \code{dnoncentchisq()} uses the (typically defining) infinite series expansion
  directly, with truncation at \code{kmax}, and terms \eqn{t_k}{t[k]} which
  are products of a Poisson probability and a central chisquare density, i.e.,
  terms \code{t.k :=  \link{dpois}(k, lambda = ncp/2) * \link{dchisq}(x, df = 2*k + df)}
  for \code{k = 0, 1, ..., kmax}.
}
\value{
  numeric vector similar to \code{x}, containing the (logged if
  \code{log=TRUE}) values of the density \eqn{f(x,*)}.
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.

  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995)
  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley.
  \cr Chapter 29, Section \emph{3  Distribution}, (29.4), p. 436.
}
\author{Martin Maechler, April 2008}
\note{
  These functions are mostly of historical interest, notably as \R's
  \code{\link{dchisq}()} was not always very accurate in the noncentral
  case, i.e., for \code{ncp > 0}.
}
\seealso{\R's own \code{\link{dchisq}()}.
}
\note{
  \R's \code{\link{dchisq}()} is typically more uniformly
  accurate than the approximations nowadays, apart from \code{dnchisqR()}
  which should behave the same.
  There may occasionally exist small differences between \code{dnchisqR(x, *)}
  and \code{\link{dchisq}(x, *)} for the same parameters.
}
\examples{
x <- sort(outer(c(1,2,5), 2^(-4:5)))
fRR <- dchisq  (x, 10, 2)
f.R <- dnchisqR(x, 10, 2)
all.equal(fRR, f.R, tol = 0) # 64bit Lnx (F 30): 1.723897e-16
stopifnot(all.equal(fRR, f.R, tol = 4e-15))
}
\keyword{math}
\keyword{distribution}
