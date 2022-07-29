\name{qnormAppr}
\title{Approximations to 'qnorm()', i.e., \eqn{z_\alpha}}
\alias{qnormAppr}
\alias{qnormUappr}
\alias{qnormUappr6}
\description{
  Relatively simple approximations to the standard normal (aka
  \dQuote{Gaussian}) quantiles, i.e., the inverse of the normal cumulative
  probability function.

  \code{qnormUappr()} is the 4-coefficient approximation to (the \bold{u}pper tail)
  standard normal quantiles, \code{\link{qnorm}()}, used in some
  \code{qbeta()} computations.

  \code{qnormUappr6()} is the \dQuote{traditional} 6-coefficient approximation to
  \code{\link{qnorm}()}, see in \sQuote{Details}.
}
\usage{% ---->>>> ../R/beta-fns.R <<<<--------
qnormAppr(p)
qnormUappr(p, lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
           lower.tail = FALSE, log.p = FALSE,
           tLarge = 1e10)
qnormUappr6(p,
            lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
               # ~= log(1-p) -- independent of lower.tail, log.p
            lower.tail=FALSE, log.p=FALSE,
            tLarge = 1e10)
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
    \eqn{\log(p)} in argument \code{p}.}
  \item{tLarge}{a large number \eqn{t0}{t_0};
    if \eqn{t >= t0}{t \ge t_0}, where \eqn{t := sqrt(-2 * lp)}{t :=
      \sqrt{-2 \cdot lp}}, the result will be \eqn{= t}.}
}
\details{
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
}
\author{Martin Maechler}
\seealso{
  \code{\link{qnorm}} (in base \R package \pkg{stats}), and importantly,
  \code{\link{qnormR}} (in our \pkg{DPQ}).
}
\examples{
pp <- c(.001, .005, .01, .05, (1:9)/10, .95, .99, .995, .999)
z_p <- qnorm(pp)
(R <- cbind(pp, z_p, qA = qnormAppr(pp),
            qUA = qnormUappr(pp, lower.tail= TRUE),
            qA6 = qnormUappr6(pp, lower.tail=TRUE)))
## Errors, absolute and relative:
cbind(pp, (relE <- cbind(
               errA  = z_p - R[,"qA" ],
               errUA = z_p - R[,"qUA"],
               rE.A  = 1 - R[,"qA" ]/z_p,
               rE.UA = 1 - R[,"qUA"]/z_p,
               rE.A6 = 1 - R[,"qA6"]/z_p)))

lp <- -c(1000, 500, 200, 100, 50, 20:10, seq(9.75, 0, by = -1/8))
qnormUappr(lp=lp) # 'p' need not be specified if 'lp' is

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
abline(h=0, v=1/2, lty=2, col="gray")

curve(qnormUappr(lp=x) / qnorm(x, log.p=TRUE) - 1,
      -2000, -.1, ylim = c(-2e-4, 1e-4), n=1001,
      main = "Relative Error of  qnormUappr(lp=x)"); mtext(" & qnormUappr6()  [log.p scale]", col=2)
curve(qnormUappr6(lp=x) / qnorm(x, log.p=TRUE) - 1, add=TRUE, col=2, n=1001)
abline(h=0, v=1/2, lty=2, col="gray")

## zoom out much more --> see jump (to be closer to truth) at ~ x = -650'000:
curve(qnormUappr(lp=x) / qnorm(x, log.p=TRUE) - 1,
      -800000, -.1, ylim = 2e-5*c(-1,1), n=1001,
      main = "Relative Error of  qnormUappr(lp=x)  [log.p scale]")
abline(h=0, v=1/2, lty=2, col="gray")
(t_ <- pretty(sqrt(-2*c(-830000, 0))))
axis(3, at=-t_^2/2, labels=t_, mgp=c(1.3, .4,0))

}
\keyword{distribution}
