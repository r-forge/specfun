\name{dnt}
\title{Non-central t-Distribution Density}
\alias{dntJKBm}
\alias{dtWVm}
\note{Package \pkg{DPQ}'s \code{\link[DPQ]{dntJKBf}()} \bold{is} already fully
  mpfr-ized, and hence \code{dntJKBm()} is redundant, and therefore deprecated.
}
\description{
  \code{dntJKBm} is a fully \CRANpkg{Rmpfr}-ified vectorized version of
  %%\link[DPQ]{dntJKBf}()} -- logically correct, but wrong by the "link to file here" rule
  \code{\link[DPQ]{dntJKBf}()} from \CRANpkg{DPQ}; see there.

  \code{dtWVm(x, df, ncp)} computes the density function \eqn{f(x)} of the t distribution with
  \code{df} degrees of freedom and non-centrality parameter \code{ncp},
  according to Wolfgang Viechtbauer's proposal in 2002, using an asymptotic
  formula for \dQuote{large} \code{df}\eqn{ = \nu}.
}
\usage{
dntJKBm(x, df, ncp, log = FALSE, M = 1000) # __ Deprecated __ use DPQ :: dntJKBf
dtWVm  (x, df, ncp, log = FALSE)
}
\arguments{
  \item{x}{numeric or \code{"\link[Rmpfr:mpfr-class]{mpfr}"} vector.}
  \item{df}{degrees of freedom (\eqn{> 0}, maybe non-integer).  \code{df = Inf} is allowed.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta};
    If omitted, use the central t distribution.}
  \item{log}{as in \code{\link{dt}()}, a logical indicating if
    \eqn{\log(f(x,*))} should be returned instead of \eqn{f(x,*)}.}
  \item{M}{the number of terms to be used, a positive integer.}
}
\details{
  See \code{\link[DPQ]{dtWV}}'s details (package \CRANpkg{DPQ}).

  As \pkg{DPQ}'s \code{\link[DPQ]{dntJKBf}()} \bold{is} already fully
  mpfr-ized, \code{dntJKBm()} is deprecated.
}
\value{
  an \code{\link[Rmpfr:mpfr-class]{mpfr}} vector of the same length as the maximum
  of the lengths of \code{x, df, ncp}.
}
\author{Martin Maechler}
\seealso{
 \R's \code{\link{dt}}, and package \CRANpkg{DPQ}'s
 \code{\link[DPQ]{dntJKBf}()} and \code{\link[DPQ]{dtWV}()}.
}
\examples{
tt <- seq(0, 10, len = 21)
ncp <- seq(0, 6, len = 31)
dt3R  <- outer(tt, ncp, dt  , df = 3)
dt3WV <- outer(tt, ncp, dtWVm, df = 3)
all.equal(dt3R, dt3WV) # rel.err 0.00063
dt25R  <- outer(tt, ncp, dt  , df = 25)
dt25WV <- outer(tt, ncp, dtWVm, df = 25)
all.equal(dt25R, dt25WV) # rel.err 1.1e-5

x <- -10:700
fx  <- dt  (x, df = 22, ncp =100)
lfx <- dt  (x, df = 22, ncp =100, log=TRUE)
lfV <- dtWVm(x, df = 22, ncp =100, log=TRUE)

head(lfx, 15) # shows that R's dt(*, log=TRUE) implementation is "quite suboptimal"

## graphics
opa <- par(no.readonly=TRUE)
par(mar=.1+c(5,4,4,3), mgp = c(2, .8,0))
plot(fx ~ x, type="l")
par(new=TRUE) ; cc <- c("red", adjustcolor("orange", 0.4))
plot(lfx ~ x, type = "o", pch=".", col=cc[1], cex=2, ann=FALSE, yaxt="n")
sfsmisc::eaxis(4, col=cc[1], col.axis=cc[1], small.args = list(col=cc[1]))
lines(x, lfV, col=cc[2], lwd=3)
dtt1 <- "      dt"; dtt2 <- "(x, df=22, ncp=100"; dttL <- paste0(dtt2,", log=TRUE)")
legend("right", c(paste0(dtt1,dtt2,")"), paste0(c(dtt1,"dtWVm"), dttL)),
       lty=1, lwd=c(1,1,3), col=c("black", cc), bty = "n")
par(opa) # reset


## For dntJKBm(), see  example(dntJKBf, package="DPQ")
}
\keyword{distribution}
\keyword{math}
