\name{qtR}
\alias{qtR}
\alias{qtR1}
\title{Pure \R Implementation of \R's C-level t-Distribution Quantiles \code{qt()}}
\description{
  A pure \R implementation of \R's Mathlib own C-level \code{\link{qt}()} function.
  \cr % in  ../R/t-dist.R
  \code{qtR()} is simply defined as \preformatted{qtR <- Vectorize(qtR1, c("p","df"))}
  where in \code{qtR1(p, df, *)} both \code{p} and \code{df} must be of length one.
}
\usage{
qtR1(p, df, lower.tail = TRUE, log.p = FALSE,
     eps = 1e-12, d1_accu = 1e-13, d1_eps = 1e-11,
     itNewt = 10L, epsNewt = 1e-14, logNewton = log.p,
     verbose = FALSE)
qtR (p, df, lower.tail = TRUE, log.p = FALSE,
     eps = 1e-12, d1_accu = 1e-13, d1_eps = 1e-11,
     itNewt = 10L, epsNewt = 1e-14, logNewton = log.p,
     verbose = FALSE)
}
\arguments{
  \item{p, df}{vectors of probabilities and degrees of freedom, see \code{\link{qt}}.}
  \item{lower.tail, log.p}{logical; see \code{\link{qt}}.}
  \item{eps}{non-negative tolerance for checking if \code{df} is
    \dQuote{very close} to \code{1} or \code{2}, respectively (when a
    special branch will be chosen).}
  \item{d1_accu, d1_eps}{non-negative tolerances only for the \code{df < 1} cases.}
  \item{itNewt}{integer, the maximal number of final Newton(-Raphson) steps.}
  \item{epsNewt}{non-negative convergence tolerance for the final Newton steps.}
  \item{logNewton}{logical, in case of \code{log.p=TRUE} indicating if
    final Newton steps should happen in log-scale.}
  \item{verbose}{logical indicating if diagnostic console output should be produced.}
}
\value{
  numeric vector of t quantiles, properly recycled in \code{(p, df)}.
}
\author{Martin Maechler}
\seealso{
 \code{\link{qtU}} and \R's \code{\link{qt}}.
}
\examples{
## Inspired from Bugzilla PR#16380
pxy <- curve(pt(-x, df = 1.09, log.p = TRUE), 4e152, 1e156, log="x", n=501)
qxy <- curve(-qt(x, df = 1.09, log.p = TRUE), -392, -385, n=501, log="y", col=4, lwd=2)
lines(x ~ y, data=pxy, col = adjustcolor(2, 1/2), lwd=5, lty=3)
## now our "pure R" version:

qRy <- -qtR(qxy$x, df = 1.09, log.p = TRUE)
all.equal(qRy, qxy$y) # "'is.NA' value mismatch: 14 in current 0 in target" for R <= 4.2.1
cbind(as.data.frame(qxy), qRy, D = qxy$y - qRy)
plot((y - qRy) ~ x, data = qxy, type="o", cex=1/4)

qtR1(.1, .1, verbose=TRUE)
pt(qtR(-390.5, df=1.10, log.p=TRUE, verbose=TRUE, itNewt = 100), df=1.10, log.p=TRUE)/-390.5 - 1
## qt(p=     -390.5, df=        1.1, *) -- general case
##  -> P=2.55861e-170, neg=TRUE, is_neg_lower=TRUE; -> final P=5.11723e-170
## usual 'df' case:  P_ok:= P_ok1 = TRUE, y=3.19063e-308, P..., !P_ok: log.p2=-390.5, y=3.19063e-308
## !P_ok && x < -36.04: q=5.87162e+153
## P_ok1: log-scale Taylor (iterated):
## it= 1, .. d{q}1=exp(lF - dt(q,df,log=T))*(lF - log(P/2)) = -5.03644e+152; n.q=5.36798e+153
## it= 2, .. d{q}1=exp(lF - dt(q,df,log=T))*(lF - log(P/2)) =  2.09548e+151; n.q=5.38893e+153
## it= 3, .. d{q}1=exp(lF - dt(q,df,log=T))*(lF - log(P/2)) =  4.09533e+148; n.q=5.38897e+153
## it= 4, .. d{q}1=exp(lF - dt(q,df,log=T))*(lF - log(P/2)) =   1.5567e+143; n.q=5.38897e+153
## [1] 0
##    === perfect!
pt(qtR(-391, df=1.10, log.p=TRUE, verbose=TRUE),
   df=1.10, log.p=TRUE)/-391 - 1 # now perfect
}
