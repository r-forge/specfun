\name{bpser}
\alias{bpser}
\title{\code{pbeta()}  'bpser' series computation}
\description{
  Compute the \code{bpser} series approximation of \code{\link{pbeta}}, the
  incomplete beta function.
  Note that when \code{b} is integer valued, the series is a \emph{sum} of
  \eqn{b+1} terms.
}
\usage{
bpser(a, b, x, log.p = FALSE, eps = 1e-15, verbose = FALSE, warn = TRUE)
}

\arguments{
  \item{a,b}{numeric and non-negative, the two shape parameters of the beta distribution.}
  \item{x}{numeric vector of abscissa values in \eqn{[0,1]}.}
  \item{log.p}{a \code{\link{logical}} if \code{\link{log}(prob)} should be
    returned, allowing to avoid underflow much farther \dQuote{out in the tails}.}
  \item{eps}{series convergence (and other) tolerance, a small positive number.}
  \item{verbose}{a \code{\link{logical}} indicating if some intermediate
    results should be printed to the console.}
  \item{warn}{a \code{\link{logical}} indicating if \code{bpser()}
    computation problems should be warned about \emph{in addition} to return
    a non-zero error code.} % ..........<<
}
%% \details{
%% }
\value{
  a \code{\link{list}} with components
 \item{r}{the resulting \code{\link{numeric}} vector.}
 \item{ier}{an integer vector of the same length as \code{x}, providing one
   error code for the computation in each \code{r[i]}.}
}
\references{
  TOMS 708, see \code{\link{pbeta}}
}
\author{Martin Maechler, ported to \pkg{DPQ}; R-Core team for the code in \R.}
\seealso{
  \R's \code{\link{pbeta}};
  \pkg{DPQ}'s \code{\link{pbetaRv1}()}, and \code{\link{Ixpq}()};
  \CRANpkg{Rmpfr}'s \code{\link[Rmpfr]{pbetaI}}
}
\examples{
with(bpser(100000, 11, (0:64)/64), # all 0 {last one "wrongly"}
     stopifnot(r == c(rep(0, 64), 1), err == 0))
bp1e5.11L <- bpser(100000, 11, (0:64)/64, log.p=TRUE)# -> 2 "underflow to -Inf" warnings!
pbe <- pbeta((0:64)/64, 100000, 11, log.p=TRUE)

## verbose=TRUE showing info on number of terms / iterations
ps11.5 <- bpser(100000, 11.5, (0:64)/64, log.p=TRUE, verbose=TRUE)
}
\keyword{distribution}

