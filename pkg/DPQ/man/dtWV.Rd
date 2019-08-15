\name{dtWV}
\alias{dtWV}
\title{Noncentral t Distribution Density}
\description{
  Compute the density function \eqn{f(x)} of the t distribution with \code{df} degrees of freedom
  and non-centrality parameter \code{ncp}.
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
  \item{log}{logical; if TRUE, \eqn{\log(f(x))} is returned instead of \eqn{f(x)}.}
}
%% \details{

%% }
\value{
  TODO
}
\references{
  Resnikoff, George J. and Lieberman, Gerald J. (1957)
  \emph{Tables of the non-central t-distribution};
  Technical report no. 32 (\code{LIE ONR 32}), April 1, 1957;
  Applied Math. and Stat. Lab., Stanford University.
  \url{https://statistics.stanford.edu/research/tables-non-central-t-distribution-density-function-cumulative-distribution-function-and}
}
\author{Wolfgang Viechtbauer (2002) post to R-help
  (\url{https://stat.ethz.ch/pipermail/r-help/2002-October/026044.html})
  as very first proposal of a noncentral t density,
  and Martin Maechler (\code{log} argument; tweaks).
}
\seealso{
 \code{\link{dt}}, \R's (C level) implementation of the (non-central) t density.
}
\examples{
tt <- seq(0, 10, len = 21)
ncp <- seq(0, 6, len = 31)
dt3R  <- outer(tt, ncp, dt  , df = 3)
dt3WV <- outer(tt, ncp, dtWV, df = 3)
all.equal(dt3R, dt3WV)

x <- -10:700
fx  <- dt(x, df = 1.1, ncp =100)
lfx <- dt  (x, df = 1.1, ncp =100, log=TRUE)
lfV <- dtWV(x, df = 1.1, ncp =100, log=TRUE)

head(lfx, 20) # shows that R's dt(*, log=TRUE) implementation is "quite suboptimal"

## graphics
opa <- par(no.readonly=TRUE)
par(mar=.1+c(5,4,4,3), mgp = c(2, .8,0))
plot(fx ~ x, type="l")
par(new=TRUE)
plot(lfx ~ x, type = "o", pch=".", col="red",cex=2, ann=FALSE, yaxt="n")
sfsmisc::eaxis(4)
lines(x, lfV, type = "o", pch=".", col=adjustcolor("orange",0.5), cex=3, lwd=3)
}
\keyword{distribution}
\keyword{math}
