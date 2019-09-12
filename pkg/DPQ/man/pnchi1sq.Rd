\name{pnchi1sq}
\title{(Probabilities of Non-Central Chisquared Distribution for Special Cases}
\alias{pnchi1sq}
\alias{pnchi3sq}%-- __end of__ >>> ../R/pnchisq.R <<<
\description{
  Computes probabilities for the non-central chi squared distribution, in
  special cases, currently for \code{df = 1} and \code{df = 3}.
}
\usage{
pnchi1sq(q, ncp = 0, lower.tail = TRUE, log.p = FALSE, epsS = .01)
pnchi3sq(q, ncp = 0, lower.tail = TRUE, log.p = FALSE, epsS = 1e-4)
}
\arguments{
  \item{q}{number ( \sQuote{quantile}, i.e., abscissa value.)}
  %\item{df}{degrees of freedom \eqn{> 0}, maybe non-integer.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{pchisq}()}.}
  \item{epsS}{small number, determining where to switch from the
    \dQuote{small case} to the regular case, namely by defining
    \code{small <- sqrt(q/ncp) <= epsS}.}
}
%% \details{
%% }
\value{
  a numeric vector \dQuote{like} \code{q+ncp}, i.e., recycled to common length.
}
\references{
  Johnson et al.(1995), see \sQuote{References} in \code{\link{pnchisqPearson}}.
}
\author{Martin Maechler, notably the Taylor approximations in the
  \dQuote{small} cases.}
\seealso{
 \code{\link{pchisq}}, the (simple and R-like) approximations, such as
 \code{\link{pnchisqPearson}} and the wienergerm approximations,
 \code{\link{pchisqW}()} etc.
}
\examples{
qq <- seq(9500, 10500, length=1000)
m1 <- cbind(pch = pchisq  (qq, df=1, ncp = 10000),
            p1  = pnchi1sq(qq,       ncp = 10000))
matplot(qq, m1, type = "l"); abline(h=0:1, v=10000+1, lty=3)
all.equal(m1[,"p1"], m1[,"pch"], tol=0) # for now,  2.37e-12

m3 <- cbind(pch = pchisq  (qq, df=3, ncp = 10000),
             p3 = pnchi3sq(qq,       ncp = 10000))
matplot(qq, m3, type = "l"); abline(h=0:1, v=10000+3, lty=3)
all.equal(m3[,"p3"], m3[,"pch"], tol=0) # for now,  1.88e-12

stopifnot(exprs = {
  all.equal(m1[,"p1"], m1[,"pch"], tol=1e-10)
  all.equal(m3[,"p3"], m3[,"pch"], tol=1e-10)
})

### Very small 'x' i.e., 'q' would lead to cancellation: -----------

##  df = 1 --------------

qS <- c(0, 2^seq(-40,4, by=1/16))
m1s <- cbind(pch = pchisq  (qS, df=1, ncp = 1)
           , p1.0= pnchi1sq(qS,       ncp = 1, epsS = 0)
           , p1.4= pnchi1sq(qS,       ncp = 1, epsS = 1e-4)
           , p1.3= pnchi1sq(qS,       ncp = 1, epsS = 1e-3)
           , p1.2= pnchi1sq(qS,       ncp = 1, epsS = 1e-2)
        )
cols <- adjustcolor(1:5, 1/2); lws <- seq(4,2, by = -1/2)
abl.leg <- function(x.leg = "topright", epsS = 10^-(4:2), legend = NULL)
{
   abline(h = .Machine$double.eps, v = epsS^2,
          lty = c(2,3,3,3), col= adjustcolor(1, 1/2))
   if(is.null(legend))
     legend <- c(quote(epsS == 0), as.expression(lapply(epsS,
                             function(K) substitute(epsS == KK,
                                                    list(KK = formatC(K, w=1))))))
   legend(x.leg, legend, lty=1:4, col=cols, lwd=lws, bty="n")
}
matplot(qS, m1s, type = "l", log="y" , col=cols, lwd=lws)
matplot(qS, m1s, type = "l", log="xy", col=cols, lwd=lws) ; abl.leg("right")
## ====  "Errors" ===================================================
## Absolute: -------------------------
matplot(qS,     m1s[,1] - m1s[,-1] , type = "l", log="x" , col=cols, lwd=lws)
matplot(qS, abs(m1s[,1] - m1s[,-1]), type = "l", log="xy", col=cols, lwd=lws)
abl.leg("bottomright")
## Relative: -------------------------
matplot(qS,     1 - m1s[,-1]/m1s[,1] , type = "l", log="x",  col=cols, lwd=lws)
abl.leg()
matplot(qS, abs(1 - m1s[,-1]/m1s[,1]), type = "l", log="xy", col=cols, lwd=lws)
abl.leg()
%% all.equal(m1s[,"p1"], m1s[,"pch"], tol=0) # for now,  2.37e-12

##  df = 3 --------------  %% FIXME:  the 'small' case is clearly wrong <<<

qS <- c(0, 2^seq(-40,4, by=1/16))
ee <- c(1e-3, 1e-2, .04)
m3s <- cbind(pch = pchisq  (qS, df=3, ncp = 1)
           , p1.0= pnchi3sq(qS,       ncp = 1, epsS = 0)
           , p1.3= pnchi3sq(qS,       ncp = 1, epsS = ee[1])
           , p1.2= pnchi3sq(qS,       ncp = 1, epsS = ee[2])
           , p1.1= pnchi3sq(qS,       ncp = 1, epsS = ee[3])
        )
matplot(qS, m3s, type = "l", log="y" , col=cols, lwd=lws)
matplot(qS, m3s, type = "l", log="xy", col=cols, lwd=lws); abl.leg("right", ee)
## ====  "Errors" ===================================================
## Absolute: -------------------------
matplot(qS,     m3s[,1] - m3s[,-1] , type = "l", log="x" , col=cols, lwd=lws)
matplot(qS, abs(m3s[,1] - m3s[,-1]), type = "l", log="xy", col=cols, lwd=lws)
abl.leg("right", ee)
## Relative: -------------------------
matplot(qS,     1 - m3s[,-1]/m3s[,1] , type = "l", log="x",  col=cols, lwd=lws)
abl.leg(, ee)
matplot(qS, abs(1 - m3s[,-1]/m3s[,1]), type = "l", log="xy", col=cols, lwd=lws)
abl.leg(, ee)
}
\keyword{distribution}
\keyword{math}
