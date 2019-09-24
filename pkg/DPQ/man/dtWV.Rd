\name{dtWV}
\alias{dtWV}
\title{Noncentral t Distribution Density by W.V.}
\description{
  Compute the density function \eqn{f(x)} of the t distribution with
  \code{df} degrees of freedom and non-centrality parameter \code{ncp},
  according to Wolfgang Viechtbauer's proposal in 2002.
}
\usage{%--> ../R/t-nonc-fn.R  + ?dt at ~/R/D/r-devel/R/src/library/stats/man/TDist.Rd
dtWV(x, df, ncp = 0, log = FALSE)
}
\arguments{
  \item{x}{numeric vector.}
  \item{df}{degrees of freedom (\eqn{> 0}, maybe non-integer).  \code{df
      = Inf} is allowed.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta};
    If omitted, use the central t distribution.}
  \item{log}{logical; if TRUE, \eqn{log(f(x))} is returned instead of \eqn{f(x)}.}
}
\details{
  The formula used is \dQuote{asymptotic}: Resnikoff and Lieberman (1957),
  p.1 and p.25ff, proposed to use recursive polynomials for (\emph{integer !})
  degrees of freedom \eqn{f = 1,2,\dots, 20}, and then, for
  \code{df}\eqn{ = f > 20}, use the asymptotic approximation which
  Wolfgang Viechtbauer proposed as a first version of a non-central t
  density for \R (when \code{\link{dt}()} did not yet have an \code{ncp}
  argument).
}
\value{
  numeric vector of density values, properly recycled in \code{(x, df, ncp)}.
}
\references{
  Resnikoff, George J. and Lieberman, Gerald J. (1957)
  \emph{Tables of the non-central t-distribution};
  Technical report no. 32 (\code{LIE ONR 32}), April 1, 1957;
  Applied Math. and Stat. Lab., Stanford University.
  \url{https://statistics.stanford.edu/research/tables-non-central-t-distribution-density-function-cumulative-distribution-function-and}
}
\author{Wolfgang Viechtbauer (2002) post to R-help
  (\url{https://stat.ethz.ch/pipermail/r-help/2002-October/026044.html}),
  and Martin Maechler (\code{log} argument; tweaks, notably recycling).
}
\seealso{
 \code{\link{dt}}, \R's (C level) implementation of the (non-central) t density.
}
\examples{
tt <- seq(0, 10, len = 21)
ncp <- seq(0, 6, len = 31)
dt3R  <- outer(tt, ncp, dt  , df = 3)
dt3WV <- outer(tt, ncp, dtWV, df = 3)
all.equal(dt3R, dt3WV) # rel.err 0.00063
dt25R  <- outer(tt, ncp, dt  , df = 25)
dt25WV <- outer(tt, ncp, dtWV, df = 25)
all.equal(dt25R, dt25WV) # rel.err 1.1e-5

x <- -10:700
fx  <- dt  (x, df = 22, ncp =100)
lfx <- dt  (x, df = 22, ncp =100, log=TRUE)
lfV <- dtWV(x, df = 22, ncp =100, log=TRUE)

head(lfx, 20) # shows that R's dt(*, log=TRUE) implementation is "quite suboptimal"

## graphics
opa <- par(no.readonly=TRUE)
par(mar=.1+c(5,4,4,3), mgp = c(2, .8,0))
plot(fx ~ x, type="l")
par(new=TRUE) ; cc <- c("red", adjustcolor("orange", 0.4))
plot(lfx ~ x, type = "o", pch=".", col=cc[1], cex=2, ann=FALSE, yaxt="n")
sfsmisc::eaxis(4, col=cc[1], col.axis=cc[1], small.args = list(col=cc[1]))
lines(x, lfV, col=cc[2], lwd=3)
dtt1 <- "      dt"; dtt2 <- "(x, df=22, ncp=100"; dttL <- paste0(dtt2,", log=TRUE)")
legend("right", c(paste0(dtt1,dtt2,")"), paste0(c(dtt1,"dtWV"), dttL)),
       lty=1, lwd=c(1,1,3), col=c("black", cc), bty = "n")
par(opa) # reset

}
\keyword{distribution}
\keyword{math}
