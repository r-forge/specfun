\name{qnormR}
\title{Pure R version of \R's \code{qnorm()} with Diagnostics and Tuning Parameters}
\alias{qnormR1}
\alias{qnormR}
\description{
  Computes \R level implementations of \R's \code{\link{qnorm}()} as
  implemented in C code (in \R's \file{Rmathlib}), historically and present.
}
\usage{
qnormR1(p, mu = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, trace = 0, version = )
qnormR (p, mu = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, trace = 0,
        version = c("4.0.x", "1.0.x", "1.0_noN", "2020-10-17", "2022-08-04"))
}
\arguments{
  \item{p}{probability \eqn{p}, \eqn{1-p}, or \eqn{\log(p)},
    \eqn{\log(1-p)}, depending on \code{lower.tail} and \code{log.p}.}
  \item{mu}{mean of the normal distribution.}
  \item{sd}{standard deviation of the normal distribution.}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qnorm}()}.}
  \item{trace}{logical or integer; if positive or \code{TRUE}, diagnostic
    output is printed to the console during the computations.}
  \item{version}{a \code{\link{character}} string specifying which version
    or variant is used.  The \emph{current} default, \code{"4.0.x"} is the
    one used in \R versions up to 4.0.x.
    The two \code{"1.0*"} versions are as used up to \R 1.0.1, based on
    Algorithm AS 111, improved by a branch for extreme tails by Wichura,
    \emph{and} a final Newton step which is only sensible when
    \code{log.p=FALSE}.  That final stepped is skipped for
    \code{version = "1.0_noN"}, \dQuote{noN} := \dQuote{no Newton}.
    %%
    \code{"2020-10-17"} is the one committed to the R development sources
    on 2020-10-17, which prevents the worst for very large \eqn{|p|} when
    \code{log.p=TRUE}.
    \code{"2022-08-04"} uses very accurate asymptotic formulas found on
    that date and provides full double precision accuracy also for extreme
    tails.}% it has been added to R's source on Feb.3, 2023 to be part of R
	   % 4.3.0 and newer
}
\details{
  For \code{qnormR1(p, ..)}, \code{p} must be of length one, whereas
  \code{qnormR(p, m, s, ..)} works vectorized in \code{p}, \code{mu}, and
  \code{sd}.  In the \pkg{DPQ} package source, \code{qnormR} is simply the result of
  \code{\link{Vectorize}(qnormR1, ...)}.
}
\value{
  a numeric vector like the input \code{q}.
}
\references{
  For \code{version}s \code{"1.0.x"} and \code{"1.0_noN"}:
  \cr
  Beasley, J.D. and Springer, S.G. (1977)
  Algorithm AS 111: The Percentage Points of the Normal Distribution.
  \emph{JRSS C (Appied Statistics)} \bold{26}, 118--121; \doi{10.2307/2346889}.

  For the asymptotic approximations used in \code{version}s newer than
  \code{"4.0.x"}, i.e., \code{"2020-10-17"} and later, see the
  reference(s) on \code{\link{qnormAsymp}}'s help page.
}
\author{Martin Maechler}
\seealso{
  \code{\link{qnorm}}, \code{\link{qnormAsymp}}.
}
\examples{
qR <- curve(qnormR, n = 2^11)
abline(h=0, v=0:1, lty=3, col=adjustcolor(1, 1/2))
with(qR, all.equal(y, qnorm(x), tol=0)) # currently shows TRUE
with(qR, all.equal(pnorm(y), x, tol=0)) # currently: mean rel. diff.: 2e-16
stopifnot(with(qR, all.equal(pnorm(y), x, tol = 1e-14)))

(ver.qn <- eval(formals(qnormR)$version)) # the possible versions
(doExtras <- DPQ:::doExtras()) # TRUE e.g. if interactive()
lp <- - 4^(1:30) # effect of  'trace = *' :
qpAll <- sapply(ver.qn, function (V)
    qnormR(lp, log.p=TRUE, trace=doExtras, version = V))
head(qpAll) # the "1.0" versions underflow quickly ..

cAdj <- adjustcolor(palette(), 1/2)
matplot(-lp, -qpAll, log="xy", type="l", lwd=3, col=cAdj, axes=FALSE,
        main = "- qnormR(lp, log.p=TRUE, version = * )")
sfsmisc::eaxis(1, nintLog=15, sub=2); sfsmisc::eaxis(2)
lines(-lp, sqrt(-2*lp), col=cAdj[ncol(qpAll)+1])
leg <- as.expression(c(paste("version=", ver.qn), quote(sqrt(-2\%*\%lp))))
matlines(-lp, -qpAll[,2:3], lwd=6, col=cAdj[2:3])
legend("top", leg, bty='n', col=cAdj, lty=1:3, lwd=2)

## Showing why/where R's qnorm() was poor up to 2020: log.p=TRUE extreme tail
##% MM: more TODO? --> ~/R/MM/NUMERICS/dpq-functions/qnorm-extreme-bad.R
qs <- 2^seq(0, 155, by=1/8)
lp <- pnorm(qs, lower.tail=FALSE, log.p=TRUE)
## The inverse of pnorm() fails BADLY for extreme tails:
## this is identical to qnorm(..) in R <= 4.0.x:
qp <- qnormR(lp, lower.tail=FALSE, log.p=TRUE, version="4.0.x")
## asymptotically correct approximation :
qpA <- sqrt(- 2* lp)
##^
col2 <- c("black", adjustcolor(2, 0.6))
col3 <- c(col2, adjustcolor(4, 0.6))
## instead of going toward infinity, it converges at  9.834030e+07 :
matplot(-lp, cbind(qs, qp, qpA), type="l", log="xy", lwd = c(1,1,3), col=col3,
        main = "Poorness of qnorm(lp, lower.tail=FALSE, log.p=TRUE)",
        ylab = "qnorm(lp, ..)", axes=FALSE)
sfsmisc::eaxis(1); sfsmisc::eaxis(2)
legend("top", c("truth", "qnorm(.) = qnormR(., \"4.0.x\")", "asymp. approx"),
       lwd=c(1,1,3), lty=1:3, col=col3, bty="n")

rM <- cbind(lp, qs, 1 - cbind(relE.qnorm=qp, relE.approx=qpA)/qs)
rM[ which(1:nrow(rM) \%\% 20 == 1) ,]
}
\keyword{distribution}
