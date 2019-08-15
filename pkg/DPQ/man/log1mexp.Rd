\name{log1mexp}
\alias{log1mexp}
\title{Compute f(a) = \eqn{\mathrm{log}}{log}(1 - \eqn{\mathrm{exp}}{exp}(-a))
  Numerically Optimally}
\description{
  Compute  f(a) = log(1 - exp(-a))   quickly and numerically accurately.
}
\usage{
log1mexp(x)
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
\keyword{math}

