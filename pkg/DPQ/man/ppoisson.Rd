\name{ppoisson}
\alias{ppoisErr}
\alias{ppoisD}
\title{Direct Computation and Errors of 'ppois()' Poisson Distribution Probabilities}
\description{
  Direct computation and errors of \code{\link{ppois}} Poisson distribution
  probabilities.
%%
%%  TODO :  say more [once "we" are sure the whole thing is worth it! ]
%%  ......
}
\usage{% --------> ../R/ppois-fn.R
ppoisD(q, lambda)
ppoisErr (lambda, ppFUN = ppoisD, iP = 1e-15, verbose = FALSE)
}
\arguments{
  \item{q}{numeric vector of non-negative integer values,
    \dQuote{quantiles} at which to evaluate \code{\link{ppois}(q, la)} and
    \code{ppFUN(q, la)}.}
  \item{lambda}{positive parameter of the Poisson distribution,
    \code{lambda}\eqn{= \lambda = E[X] = Var[X]} where
    \eqn{X \sim Pois(\lambda)}{X ~ Pois(\lambda)}.}
  \item{ppFUN}{alternative \code{\link{ppois}} evaluation, by default the
    \bold{d}irect summation of \code{\link{dpois}(k, lambda)}.}
  \item{iP}{small number, \code{iP << 1}, used to construct the abscissa values
    \code{x} at which to evaluate and compare \code{ppois()} and \code{ppFUN()}:
    \code{x <- 0:qpois(1-iP, lambda = lambda)}.}
  \item{verbose}{\code{\link{logical}} indicating if extra information
    should be printed.}
}
%% \details{
%%   ....
%% }
\value{
  \code{ppoisD()} contains the poisson probabilities along \code{q}, i.e.,
  is a numeric vector of length \code{length(q)}.

  \code{re <- ppoisErr()} returns the relative \dQuote{error} of
  \code{\link{ppois}(x0, lambda)} where \code{ppFUN(x0, lambda)} is
  assumed to be the truth and \code{x0} the \dQuote{worst case}, i.e., the
  value (among \code{x}) with the largest such difference.

  Additionally, \code{attr(re, "x0")} contains that value \code{x0}.
}
\author{Martin Maechler, March 2004; 2019 ff}
\seealso{
  \code{\link{ppois}}
}
\examples{
(lams <- outer(c(1,2,5), 10^(0:3)))# 10^4 is already slow!
system.time(e1 <- sapply(lams, ppoisErr))
e1 / .Machine$double.eps

## Try another 'ppFUN' :---------------------------------
## this relies on the fact that it's *only* used on an 'x' of the form  0:M :
ppD0 <- function(x, lambda) cumsum(dpois(x, lambda=lambda))
## and test it:
p0 <- ppD0 (0:1000, lambda=10)
p1 <- ppois(0:1000, lambda=10)
stopifnot(all.equal(p0,p1, tol=8*.Machine$double.eps))

## Note that it is "infinitely" much faster but possibly less accurate
##  ("possibly": need a {pkg 'Rmpfr'}-enabled version to get that) :
system.time(p15     <- ppoisErr(2^15))               ## 3.8.. sec !
system.time(p15.0.  <- ppoisErr(2^15, ppFUN = ppD0)) ## 0.009
c(p15, p15.0.) / .Machine$double.eps # on Lnx 64b, see  (-10  2.5)

## lapply(), so you see "x0" values :
str(e0. <- lapply(lams, ppoisErr, ppFUN = ppD0))

## The first version [called 'err.lambd0()' for years] used simple  cumsum(dpois(..))
## NOTE: It is *stil* much faster, as it relies on special  x == 0:M  relation
## Author: Martin Maechler, Date:  1 Mar 2004, 17:40
##
## Recreate that version as  'ppoisErr0()'
## (and learn about smart use of `formals<-` and `environment<-`) :
ppoisErr0 <- ppoisErr
formals(ppoisErr0)[["ppFUN"]] <- function(x, lambda) cumsum(dpois(x, lambda=lambda))
environment(ppoisErr0) <- asNamespace("DPQ")
e0 <- sapply(lams, ppoisErr0)
stopifnot(e0 == unlist(e0.))
all.equal(e1, e0) # typically TRUE,  though small "random" differences:
cbind(e1, e0) * 2^53 # on Lnx 64b, seeing integer values in {-24, .., 33}
}
\keyword{distribution}
