\name{numer-utils}
\title{Numerical Utilities - Functions, Constants}
%% Numeric Constants .. ---------------------------
\alias{M_LN2}
\alias{M_SQRT2}
\alias{M_minExp}
\alias{M_cutoff}
\alias{G_half}
%% Functions ... ----------------------------------
\alias{all_mpfr}
\alias{any_mpfr}
\alias{logr}
\alias{modf}
\alias{okLongDouble}
%--------------------------------------------------
\description{
  The \pkg{DPQ} package provides some numeric constants used in some of its
  distribution computations.

  \code{all_mpfr()} and \code{any_mpfr()} return \code{\link{TRUE}}
  iff all (or \sQuote{any}, respectively) of their arguments inherit from
  \code{\link{class}} \code{"mpfr"} (from package \CRANpkg{Rmpfr}).

  \code{logr(x,a)} computes \code{\link{log}(x / (x + a))} in a numerically
  stable way.

  \code{modf(x)} splits each \code{x} into integer part (as
  \code{\link{trunc}(x)}) and fractional (remainder) part in \eqn{(-1, 1)}
  and corresponds to the \R version of the C99 (and POSIX) standard C (and C++) mathlib
  functions of the same name.
}
\usage{
## Numeric Constants : % mostly in   ../R/beta-fns.R
M_LN2        # = log(2)  = 0.693....
M_SQRT2      # = sqrt(2) = 1.4142...
M_cutoff     # := If |x| > |k| * M_cutoff, then  log[ exp(-x) * k^x ]  =~=  -x
             #  = 3196577161300663808 ~= 3.2e+18
M_minExp     # = log(2) * .Machine$double.min.exp # ~= -708.396..
G_half       # = sqrt(pi) = Gamma( 1/2 )

## Functions :
all_mpfr(...)
any_mpfr(...)
logr(x, a)    # == log(x / (x + a)) -- but numerically smart; x >= 0, a > -x
modf(x)
okLongDouble(lambda = 999, verbose = 0L, tol = 1e-15)
}
\arguments{
  \item{\dots}{numeric or \code{"mpfr"} numeric vectors.}
  \item{x, a}{number-like, not negative, now may be vectors of
    \code{length(.) > 1}.}
  \item{lambda}{a number, typically in the order of 500--10'000.}
  \item{verbose}{a non-negative integer, if not zero,
    \code{okLongDouble()} prints the intermediate long double
    computations' results.}
  \item{tol}{numerical tolerance used to determine the accuracy required
    for near equality in \code{okLongDouble()}.}
}
\value{
  The numeric constant in the first case; a numeric (or "mpfr") vector of
  appropriate size in the 2nd case.

  \code{okLongDouble()} returns a \code{\link{logical}},
  \code{\link{TRUE}} iff the long double arithmetic with \code{expl()} and
  \code{logl()} seems to work accurately
  and consistently for \code{exp(-lambda)} and \code{log(lambda)}.
}
\details{
  \describe{
    \item{\code{all_mpfr()},}{}
    \item{\code{all_mpfr()} :}{test if \code{\link{all}} or
      \code{\link{any}} of their arguments or of class \code{"mpfr"} (from
      package \CRANpkg{Rmpfr}).  The arguments are evaluated only until
      the result is determined, see the example.}

    \item{\code{logr()}}{computes \eqn{\log( x / (x+a) )}{log(x / (x+a))} in a numerically
      stable way.}%% which means .... (TODO)

    %% \item{\code{..()}}{ ..}
    %% \item{\code{..()}}{ ..}
  }
}
%% \references{
%%   ...
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{.Machine}}
}
\examples{
(Ms <- ls("package:DPQ", pattern = "^M"))
lapply(Ms, function(nm) { cat(nm,": "); print(get(nm)) }) -> .tmp

logr(1:3, a=1e-10)

okLongDouble() # typically TRUE, but not e.g. in a valgrinded R-devel of Oct.2019
## Here is typically the "boundary":
rr <- uniroot(function(x) okLongDouble(x) - 1/2, c(11350, 11400), tol=1e-7)
str(rr, digits=9) ## seems somewhat platform dependent: now see
## $ root      : num 11376.563
## $ estim.prec: num 9.313e-08
## $ iter      : int 29

set.seed(2021); x <- runif(100, -7,7) 
mx <- modf(x)
with(mx, head( cbind(x, i=mx$i, fr=mx$fr) )) # showing the first cases
with(mx, stopifnot(   x == fr + i,
                      i == trunc(x),
               sign(fr) == sign(x)))
}
\keyword{math}
