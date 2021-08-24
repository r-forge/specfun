\name{log1mexp}
\title{Compute \eqn{\mathrm{log}}{log}(1 - \eqn{\mathrm{exp}}{exp}(-a)) and
  \eqn{\log(1 + \exp(x))}{log(1 + exp(x))}   Numerically Optimally}
\alias{log1mexp}
\alias{log1mexpC}
\alias{log1pexpC}
\description{
  Compute  f(a) = log(1 - exp(-a))   quickly and numerically accurately.

  \code{log1mexp()} is simple pure \R code;\cr
  \code{log1mexpC()} is an interface to \R C API (\file{Mathlib} / \file{Rmath.h})
  function.

  \code{log1pexpC()} is an interface to \R's \file{Mathlib} \code{double}
  function \code{log1pexp()} which computes \eqn{\log(1 + \exp(x))}{log(1 + exp(x))},
  accurately, notably for large \eqn{x}, say, \eqn{x > 720}.
}
\usage{
log1mexp (x)
log1mexpC(x)
log1pexpC(x)
}
\arguments{
  \item{x}{numeric vector of positive values.}
}
\author{Martin Maechler}
\seealso{% ~/R/Pkgs/copula/man/log1mexp.Rd
  The \code{\link{log1mexp}()} function in CRAN package \CRANpkg{copula},
  and the corresponding vignette (in the \sQuote{References}).
}
\references{%  ~/R/Pkgs/Rmpfr/vignettes/log1mexp-note.Rnw
  Martin \enc{Mächler}{Maechler} (2012).
  Accurately Computing \eqn{\log(1-\exp(-|a|))};
  \url{https://CRAN.R-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}.
  %  see also <<log1mexp>> in ../inst/doc/Frank-Rmpfr.Rnw
}
\examples{
l1m.xy <- curve(log1mexp(x), -10, 10, n=1001)
stopifnot(with(l1m.xy, all.equal(y, log1mexpC(x))))

x <- seq(0, 710, length=1+710*2^4); stopifnot(diff(x) == 1/2^4)
l1pm <- cbind(log1p(exp(x)),
              log1pexpC(x))
matplot(x, l1pm, type="l", log="xy") # both look the same
iF <- is.finite(l1pm[,1])
stopifnot(all.equal(l1pm[iF,2], l1pm[iF,1], tol=1e-15))
}
\keyword{math}
