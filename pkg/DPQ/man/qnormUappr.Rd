\name{qnormAppr}
\alias{qnormAppr}
\alias{qnormUappr}
\title{Approximations to 'qnorm()', i.e., \eqn{z_\alpha}}
\description{
  Relatively simple approximations to the standard normal (aka
  \dQuote{Gaussian}) quantiles, i.e., the inverse of the normal cumulative
  probability function.

  \code{qnormUappr()} is a simple approximation to (the \bold{u}pper tail)
  standard normal quantiles, \code{\link{qnorm}()}.
}
\usage{
qnormAppr(p)
qnormUappr(p, lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
           lower.tail = FALSE, log.p = FALSE)
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
}
\details{
  \code{qnormAppr(p)} uses the simple 4 coefficient rational approximation
  to \code{\link{qnorm}(p)}, % provided by ... literature ?
  to be used \emph{only} for \eqn{p > 1/2} in \code{\link{qbeta}()} computations, e.g.,
  \code{\link{qbeta.R}}.
  \cr
  The relative error of this approximation is quite \emph{asymmetric}: It
  is mainly < 0.

  \code{qnormUappr(p)} uses the same rational approximation directly for the
  \bold{U}pper tail where it is relatively good, and for the lower tail via
  \dQuote{swapping the tails}, so it is good there as well.
}
\value{
  numeric vector of (approximate) normal quantiles corresponding to
  probabilities \code{p}
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
  \code{\link{qnorm}}.
}
\examples{
pp <- c(.001, .005, .01, .05, (1:9)/10, .95, .99, .995, .999)
z_p <- qnorm(pp)
(R <- cbind(pp, z_p, qA = qnormAppr(pp), qUA = qnormUappr(pp, lower.tail=TRUE)))
## Errors, absolute and relative:
cbind(pp, (relE <- cbind(
               errA  = z_p - R[,"qA" ],
               errUA = z_p - R[,"qUA"],
               rE.A  = 1 - R[,"qA" ]/z_p,
               rE.UA = 1 - R[,"qUA"]/z_p)))

lp <- -c(1000, 500, 200, 100, 50, 20:10, seq(9.75, 0, by = -1/8))
qnormUappr(lp=lp) # 'p' need not be specified if 'lp' is

curve(qnorm(x, lower.tail=FALSE), n=1001)
curve(qnormUappr(x), add=TRUE,    n=1001, col = adjustcolor("red", 1/2))

curve(qnorm(x, lower.tail=FALSE) - qnormUappr(x), n=1001)
}
\keyword{distribution}
