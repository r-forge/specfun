\name{DPQmpfr-utils}
\alias{ldexp}
\title{Numeric / Mpfr Utilities for DPQmpfr}
\description{
  Utilities for package \CRANpkg{DPQmpfr}
}
\usage{
ldexp(f, E)
}
\arguments{% see ../../DPQ/man/fr_ld_exp.Rd
  \item{f}{\sQuote{fraction}, as such with absolute value in \eqn{[0.5, 1)}, but can
    be any numbers.}
  \item{E}{integer-valued exponent(s).}
}
\details{
  \code{ldexp()} is a simple wrapper, either calling
  \code{DPQ::\link[DPQ]{ldexp}} from \CRANpkg{DPQ} or
  \code{\link[Rmpfr]{ldexpMpfr}} from the \CRANpkg{Rmpfr} package,
  \deqn{ ldexp(f, E) := f \times 2^E,}{
         ldexp(f, E) := f * 2^E,}
  computed accurately and fast on typical platforms with internally binary
  arithmetic.
}
\value{
  % \code{ldexp()}:
  either a numeric or a \code{"\link{mpfr}"}, depending on the type of
  \code{f}, vector as (the recyled) combination of \code{f} and \code{E}.
}
\seealso{
  \code{\link[DPQ]{ldexp}} from package \CRANpkg{DPQ} and
  \code{\link[Rmpfr]{ldexpMpfr}} from package \CRANpkg{Rmpfr}.
}
\examples{
ldexp(1:10, 2)
ldexp(Rmpfr::Const("pi", 96), -2:2) # =  pi * (1/4  1/2  1  2  4)
}
\keyword{arithmetic}
