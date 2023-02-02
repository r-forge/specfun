\name{qnormAppr}
\title{Approximations to 'qnorm()', i.e., \eqn{z_\alpha}}
\alias{qnormAppr}% <- deprecated!
\alias{qnormUappr}
\alias{qnormUappr6}
\alias{qnormCappr}
\description{
  Approximations to the standard normal (aka \dQuote{Gaussian}) quantiles,
  i.e., the inverse of the normal cumulative probability function.

  The \code{qnormUappr*()} are relatively simple approximations from
  Abramowitz and Stegun, computed by Hastings(1955):
  \code{qnormUappr()} is the 4-coefficient approximation to (the \bold{u}pper tail)
  standard normal quantiles, \code{\link{qnorm}()}, used in some
  \code{qbeta()} computations.

  \code{qnormUappr6()} is the \dQuote{traditional} 6-coefficient approximation to
  \code{\link{qnorm}()}, see in \sQuote{Details}.
}
\usage{% ---->>>> ../R/beta-fns.R <<<<--------
qnormUappr(p, lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
           lower.tail = FALSE, log.p = missing(p),
           tLarge = 1e10)
qnormUappr6(p,
            lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
               # ~= log(1-p) -- independent of lower.tail, log.p
            lower.tail = FALSE, log.p = missing(p),
            tLarge = 1e10)

qnormCappr(p, k = 1) ## *implicit* lower.tail=TRUE, log.p=FALSE  >>> TODO: add! <<

qnormAppr(p) # << deprecated; use qnormUappr(..) instead!
}
\arguments{
  \item{p}{numeric vector of probabilities, possibly transformed, depending
    on \code{log.p}.  Does not need to be specified, if \code{lp} is instead.}
  \item{lp}{\code{log(1 - p*)}, assuming \eqn{p*} is the \code{lower.tail=TRUE,
      log.p=FALSE} version of \code{p}.  If passed as argument, it can be
    much more accurate than when computed from \code{p} by default.}
  \item{lower.tail}{logical; if TRUE (\emph{not} the default here!), probabilities are
    \eqn{P[X \le x]}, otherwise (by default) upper tail probabilities, \eqn{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities \eqn{p} are given as
    \eqn{\log(p)} in argument \code{p}.  Note that it is \emph{not used},
    when \code{missing(p)} and \code{lp} is specified.}
  \item{tLarge}{a large number \eqn{t0}{t_0};
    if \eqn{t >= t0}{t \ge t_0}, where \eqn{t := sqrt(-2 * lp)}{t :=
      \sqrt{-2 \cdot lp}}, the result will be \eqn{= t}.}
  \item{k}{positive integer, specifying the iterative plugin \sQuote{order}.}
}
\details{
  This is now \emph{deprecated}; use  \code{qnormUappr()} instead!
  \code{qnormAppr(p)} uses the simple 4 coefficient rational approximation
  to \code{\link{qnorm}(p)}, provided by Abramowitz and Stegun (26.2.22), p.933,
  to be used \emph{only} for \eqn{p > 1/2} and typically
  \code{\link{qbeta}()} computations, e.g., \code{\link{qbeta.R}}.
  \cr
  The relative error of this approximation is quite \emph{asymmetric}: It
  is mainly < 0.

  \code{qnormUappr(p)} uses the same rational approximation directly for the
  \bold{U}pper tail where it is relatively good, and for the lower tail via
  \dQuote{swapping the tails}, so it is good there as well.

  \code{qnormUappr6(p, *)} uses the 6 coefficient rational approximation
  to \code{\link{qnorm}(p, *)}, from Abramowitz and Stegun (26.2.23), again
  mostly useful in the outer tails.

  \code{qnormCappr(p, k)} inverts formula (26.2.24) of Abramowitz and Stegun,
  and for \eqn{k \ge 2}{k >= 2} improves it, by iterative recursive
  plug-in, using A.&S. (26.2.25).
}
\value{
  numeric vector of (approximate) normal quantiles corresponding to
  probabilities \code{p}
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.

  Hastings jr., Cecil (1955)
  \emph{Approximations for Digital Computers}.
  Princeton Univ. Press.
}
\author{Martin Maechler}
\seealso{
  \code{\link{qnorm}} (in base \R package \pkg{stats}), and importantly,
  \code{\link{qnormR}} and \code{\link{qnormAsymp}()} in this package (\pkg{DPQ}).
}
\examples{
pp <- c(.001, .005, .01, .05, (1:9)/10, .95, .99, .995, .999)
z_p <- qnorm(pp)
assertDeprecation <- function(expr, verbose=TRUE)
  tools::assertCondition(expr, verbose=verbose, "deprecatedWarning")
assertDeprecation(qA <- qnormAppr(pp))
(R <- cbind(pp, z_p, qA,
            qUA = qnormUappr(pp, lower.tail= TRUE),
            qA6 = qnormUappr6(pp, lower.tail=TRUE)))
## Errors, absolute and relative:
relEr <- function(targ, curr) { ## simplistic "smart" rel.error
    E <- curr - targ
    r <- E/targ  # simple, but fix 0/0:
    r[targ == 0 & E == 0] <- 0
    r
}
mER <- cbind(pp,
             errA  = z_p - R[,"qA" ],
             errUA = z_p - R[,"qUA"],
             rE.A  = relEr(z_p, R[,"qA" ]),
             rE.UA = relEr(z_p, R[,"qUA"]),
             rE.A6 = relEr(z_p, R[,"qA6"]))
signif(mER)

lp <- -c(1000, 500, 200, 100, 50, 20:10, seq(9.75, 0, by = -1/8))
signif(digits=5, cbind(lp # 'p' need not be specified if 'lp' is !
    , p.  = -expm1(lp)
    , qnU = qnormUappr (lp=lp)
    , qnU6= qnormUappr6(lp=lp)
    , qnA1= qnormAsymp(lp=lp, lower.tail=FALSE, order=1)
    , qnA5= qnormAsymp(lp=lp, lower.tail=FALSE, order=5)
    , qn  = qnorm(lp, log.p=TRUE)
      )) ## oops! shows *BUG* for last values where qnorm() > 0 !

curve(qnorm(x, lower.tail=FALSE), n=1001)
curve(qnormUappr(x), add=TRUE,    n=1001, col = adjustcolor("red", 1/2))

## Error curve:
curve(qnormUappr(x) - qnorm(x, lower.tail=FALSE), n=1001,
      main = "Absolute Error of  qnormUappr(x)")
abline(h=0, v=1/2, lty=2, col="gray")

curve(qnormUappr(x) / qnorm(x, lower.tail=FALSE) - 1, n=1001,
      main = "Relative Error of  qnormUappr(x)")
 abline(h=0, v=1/2, lty=2, col="gray")

curve(qnormUappr(lp=x) / qnorm(x, log.p=TRUE) - 1, -200, -1, n=1001,
      main = "Relative Error of  qnormUappr(lp=x)"); mtext(" & qnormUappr6()  [log.p scale]", col=2)
curve(qnormUappr6(lp=x) / qnorm(x, log.p=TRUE) - 1, add=TRUE, col=2, n=1001)
abline(h=0, lty=2, col="gray")

curve(qnormUappr(lp=x) / qnorm(x, log.p=TRUE) - 1,
      -2000, -.1, ylim = c(-2e-4, 1e-4), n=1001,
      main = "Relative Error of  qnormUappr(lp=x)"); mtext(" & qnormUappr6()  [log.p scale]", col=2)
curve(qnormUappr6(lp=x) / qnorm(x, log.p=TRUE) - 1, add=TRUE, col=2, n=1001)
abline(h=0, lty=2, col="gray")

## zoom out much more - switch x-axis {use '-x'} and log-scale:
curve(qnormUappr6(lp=-x) / qnorm(-x, log.p=TRUE) - 1,
      .1, 1.1e10, log = "x", ylim = 2.2e-4*c(-2,1), n=2048,
      main = "Relative Error of  qnormUappr6(lp = -x)  [log.p scale]") -> xy.q
abline(h=0, lty=2, col="gray")

## 2023-02: qnormUappr6() can be complemented with
## an approximation around center p=1/2: qnormCappr()
p <- seq(0,1, by=2^-10)
M <- cbind(p, qn=(qn <- qnorm(p)),
           reC1 = relEr(qn, qnormCappr(p)),
           reC2 = relEr(qn, qnormCappr(p, k=2)),
           reC3 = relEr(qn, qnormCappr(p, k=3)),
           reU6 = relEr(qn, qnormUappr6(p,lower.tail=TRUE)))
matplot(M[,"p"], M[,-(1:2)], type="l", col=2:7, lty=1, lwd=2,
        ylim = c(-.004, +1e-4), xlab=quote(p), ylab = "relErr")
abline(h=0, col="gray", lty=2)
oo <- options(width=99)
summary(    M[,-(1:2)])
summary(abs(M[,-(1:2)]))
options(oo)
}
\keyword{distribution}
