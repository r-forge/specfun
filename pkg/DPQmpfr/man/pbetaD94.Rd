\name{betaD94}
\title{Ding(1994) (non-central) Beta Distribution Functions}
\alias{dbetaD94}
\alias{pbetaD94}
\alias{qbetaD94}
\description{% ../R/pbeta-Ding94.R
  The three functions "p" (cumulative distribution, CDF), "d" (density
  (PDF)), and "q" (quantile) use Ding(1994)'s algorithm A, B, and C, respectively,
  each of which implements a recursion formula using only simple
  arithmetic and \code{\link{log}} and \code{\link{exp}}.

  These are particularly useful also for using with high precision
  \code{"mpfr"} numbers from the \CRANpkg{Rmpfr} CRAN package.
}
\usage{
dbetaD94(x, shape1, shape2, ncp = 0, log = FALSE,
         eps = 1e-10, itrmax = 100000L, verbose = FALSE)
pbetaD94(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
	 log_scale = (a * b > 0) && (a + b > 100 || c >= 500),
         eps = 1e-10, itrmax = 100000L, verbose = FALSE)

qbetaD94(p, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
	 log_scale = (a * b > 0) && (a + b > 100 || c >= 500),
         delta = 1e-6,
         eps = delta^2,
         itrmax = 100000L,
         iterN = 1000L,
         verbose = FALSE)
}

\arguments{
  \item{x,q}{numeric vector of values in \eqn{[0,1]} as beta variates.}
  \item{shape1, shape2}{the two shape parameters of the beta
    distribution, must be positive.}
  \item{ncp}{the noncentrality parameter; by default zero for the
    (\emph{central}) beta distribution; if positive, we have a noncentral
    beta distribution.}
  \item{p}{numeric vector of probabilities, \code{\link{log}()}ged in
    case \code{log.p} is true.}
  \item{log, log.p}{logical indicating if the density or probability
    values should be \code{\link{log}()}ged.}
  \item{lower.tail}{logical indicating if the lower or upper tail
    probability should be computed, or for \code{qbeta*()} are provided.}
  \item{eps}{a non-negative number specifying the desired accuracy for
    computing F() and f().}
  \item{itrmax}{the maximal number of steps for computing F() and f().}
  \item{delta}{[For \code{qbeta*()}:] non-negative number indicating the
    desired accuracy for computing \eqn{x_p} (the root of \eqn{pbeta*() == p}),
    i.e., the convergence tolerance for the Newton iterations.   This sets
    default \code{eps = delta^2} which is sensible but may be too small,
    such that \code{eps} should be specified in addition to \code{delta}.}
  \item{iterN}{[For \code{qbeta*()}:] The maximal number of \bold{N}ewton iterations.}
  \item{log_scale}{logical indicating if most of the computations should
    happen in \code{\link{log}} scale, which protects from \dQuote{early}
    overflow and underflow but takes more computations.  The current
    default is somewhat \emph{arbitrary}, still derived from the facts that
    \code{\link{gamma}(172)} overflows to \code{Inf} already and
    \code{\link{exp}(-750)} underflows to \code{0} already.}
  \item{verbose}{logical (or integer) indicating the amount of
    diagnostic output during computation; by default none.}
}
%% \details{
%%
%% }
\value{
  In all three cases, a numeric vector with the same attributes as
  \code{x} (or \code{q} respectively),
  containing (an approximation) to the correponding beta distribution function.
}
\references{
  Cherng G. Ding (1994) % ~/save/papers/Numerics/Ding_1994_noncentr-beta.pdf
  On the computation of the noncentral beta distribution.
  \emph{Computational Statistics & Data Analysis} \bold{18}, 449--455.
}
\author{Martin Maechler, notably \code{log_scale} was not part of Ding's
  proposals.}
\seealso{
  \code{\link{pbeta}}.  Package \CRANpkg{Rmpfr}'s \code{\link[Rmpfr]{pbetaI}()} needs
  both \code{shape1} and \code{shape2} to be integer but is typically more
  efficient than the current \code{pbetaD94()} implementation.
}
\examples{%% see also ../tests/beta-Ding94.R
          %%          ~~~~~~~~~~~~~~~~~~~~~~
## Low precision (eps, delta) values as "e.g." in Ding(94): ------------------

## Compare with  Table 3  of  Baharev_et_al 2017 %% ===> ./qbBaha2017.Rd <<<<<<<<<<<
aa <- c(0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 25)
bb <- c(1:15, 10*c(2:5, 10, 25, 50))
utime <-
 qbet <- matrix(NA_real_, length(aa), length(bb),
                dimnames = list(a = formatC(aa), b = formatC(bb)))
(doExtras <- DPQmpfr:::doExtras())
if(doExtras) qbetL <- utimeL <- utime

p <- 0.95
delta <- 1e-4
eps   <- 1e-6
system.t.usr <- function(expr)
  system.time(gcFirst = FALSE, expr)[["user.self"]]

system.time(
for(ia in seq_along(aa)) {
    a <- aa[ia]; cat("\n--==--\na=",a,":\n")
    for(ib in seq_along(bb)) {
        b <- bb[ib]; cat("\n>> b=",b,"\n")
        utime [ia, ib] <- system.t.usr(
          qbet[ia, ib] <-   qbetaD94(p, a, b, ncp = 0, delta=delta, eps=eps, verbose = 2))
        if(doExtras)
          utimeL[ia, ib] <- system.t.usr(
           qbetL[ia, ib] <-   qbetaD94(p, a, b, ncp = 0, delta=delta, eps=eps,
                                       verbose = 2, log_scale=TRUE))
    }
    cat("\n")
}
)# system.time(.): ~ 1 sec (lynne i7-7700T, Fedora 32, 2020)
sum(print(table(round(1000*utime)))) # lynne .. :
##  0  1  2  3  4  5  6  7  8  9 10 11 14 15 16 29
## 53 94 15  3  3 12  2  2  2  2  1  2  3  1  2  1
## [1] 198
if(doExtras) print(sum(print(table(round(1000*utimeL))))) # lynne .. :
}
\keyword{distribution}
\keyword{math}
