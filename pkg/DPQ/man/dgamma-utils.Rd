\name{dgamma-utils}
\title{Utility Functions for \code{dgamma()} -- Pure R Versions}
\alias{bd0}
\alias{dpois_raw}
\alias{stirlerr}
\description{
  R transcriptions of the C code utility functions for
  \code{\link{dgamma}()} and similar \dQuote{base} density functions by
  Catherine Loader.
}
\usage{
dpois_raw(x, lambda, log)
bd0(x, np, verbose = getOption("verbose"))
stirlerr(n)
}
\arguments{
  \item{x, n}{a number (\code{\link{numeric}}).}
  \item{lambda, np}{a number (\code{\link{numeric}}; distrubution parameters.}
  \item{log}{logical indicating if the log-density should be returned,
    otherwise the density at \code{x}.}
  \item{verbose}{logical indicating if some information about the
    computations are to be printed.}
}
\value{
a number
}
%% \details{
%% }
%% \references{
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{dgamma}}, % ~/R/D/r-devel/R/src/library/stats/man/GammaDist.Rd
  \code{\link{dpois}}.
}
\examples{
n <- seq(1, 50, by=1/4)
plot(n, sapply(n, stirlerr), type = "b", log="xy")
}
\keyword{distribution}

