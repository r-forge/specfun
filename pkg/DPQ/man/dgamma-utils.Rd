\name{dgamma-utils}
\title{Utility Functions for \code{dgamma()} -- Pure R Versions}
\alias{bd0}
\alias{ebd0}
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
ebd0(x, M, verbose = getOption("verbose"))
stirlerr(n, scheme = c("R3", "R4.1"),
         cutoffs = switch(scheme
                        , R3   = c(15, 35, 80, 500)
                        , R4.1 = c(15, 40, 85, 600) ## FIXME - unfinished
                          ))
}
\arguments{
  \item{x, n}{a number (\code{\link{numeric}}).}
  \item{lambda, np, M}{a number (\code{\link{numeric}}; distrubution parameters.}
  \item{log}{logical indicating if the log-density should be returned,
    otherwise the density at \code{x}.}
  \item{verbose}{logical indicating if some information about the
    computations are to be printed.}
  \item{scheme}{a \code{\link{character}} string specifying the
    \code{cutoffs} scheme.}
  \item{cutoffs}{a (increasing) numeric vector, currently required to start
    with \code{15} specifying the cutoffs to switch from 2 to 3 to 4 and 5
    terms approximations for non-small \code{n}, where the direct formula
    loses precision.}
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
st.n <- stirlerr(n) # now vectorized
stopifnot(identical(st.n, sapply(n, stirlerr)))
plot(n, st.n, type = "b", log="xy", ylab = "stirlerr(n)")

x <- 800:1200
bd0x1k <- bd0(x, np = 1000)
plot(x, bd0x1k, type="l", ylab = "bd0(x, np=1000)")
if(FALSE) ## still fails --- ebd0() is really unfinished !!
ebd0x1k <- sapply(x, function(.) ebd0(., 1000)) ## not yet vectorized
}
\keyword{distribution}

