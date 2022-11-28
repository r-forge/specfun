\name{qnormAsymp}
\alias{qnormAsymp}
\title{Asymptotic Approximation to Outer Tail of qnorm()}
\description{
  Implementing new asymptotic tail approximations of normal quantiles,
  i.e., the \R function \code{\link{qnorm}()}, mostly useful when
  \code{log.p=TRUE} and log-scale \code{p} is relatively large negative,
  i.e., \eqn{p \ll -1}{p << -1}.
}
\usage{% ---->>>> ../R/norm_f.R <<<<--------
qnormAsymp(p,
           lp = .DT_Clog(p, lower.tail = lower.tail, log.p = log.p),
           order, lower.tail = TRUE, log.p = missing(p))
}
\arguments{
  \item{p}{numeric vector of probabilities, possibly transformed, depending
    on \code{log.p}.  Does not need to be specified, if \code{lp} is instead.}
  \item{lp}{numeric (vector) of \code{log(1-p)} values; if not specified,
    computed from \code{p}, depending on \code{lower.tail} and \code{log.p}.}
  \item{order}{an integer in \eqn{\{0,1,\dots,5\}}{{0,1,..,5}}, specifying the
    approximation order.}
  \item{lower.tail}{logical; if true, probabilities are \eqn{P[X \le x]},
    otherwise upper tail probabilities, \eqn{P[X > x]}.}
  \item{log.p}{logical; if \code{TRUE} (as typical here!), probabilities
    \eqn{p} are given as \eqn{\log(p)} in argument \code{p}.}
}
\details{
  These \emph{asymptotic} approximations have been derived by Maechler(2022)
  via iterative plug-in to the well known asymptotic approximations of
  \eqn{Q(x) = 1 - \Phi(x)} from Abramowitz and Stegun (26.2.13), p.932,
  which are provided in our package \pkg{DPQ} as \code{\link{pnormAsymp}()}.
  %%
  They will be used in R >= 4.3.0's \code{qnorm()} to provide very accurate
  quantiles in the extreme tails.
}
\value{
  a numeric vector like \code{p} or \code{lp} if that was specified instead.

  The simplemost (for extreme tails) is \code{order = 0}, where the
  asymptotic approximation is simply \eqn{\sqrt{-2s}}{sqrt(-2s)} and
  \eqn{s} is \code{-lp}.
}
\references{
  Martin Maechler (2022). Asymptotic Tail Formulas For Gaussian Quantiles for R;
  \CRANpkg{DPQ} vignette, see \url{https://CRAN.R-project.org/package=DPQ}.
}
\author{Martin Maechler}
\seealso{
  The upper tail approximations in Abramowitz & Stegun, in \pkg{DPQ}
  available as \code{qnormUappr()} and \code{\link{qnormUappr6}()},
  are less accurate than our \code{order >= 1} formulas in the tails.
}
\examples{%% see also ./qnormUappr.Rd and notably ../tests/pqnorm_extreme.R <<<
lp <- -c(head(c(outer(c(5,2,1), 10^(18:1))), -2), 20:10, seq(9.75, 2, by = -1/8))
qnU6 <- qnormUappr6(lp=lp) # 'p' need not be specified if 'lp' is
qnAsy <- sapply(0:5, function(ord) qnormAsymp(lp=lp, lower.tail=FALSE, order=ord))
matplot(-lp, cbind(qnU6, qnAsy), type = "b", log = "x", pch=1:7)# all "the same"
legend("center", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:6, lty=1:5, pch=1:7) # as in matplot()

p.ver <- function() mtext(R.version.string, cex=3/4, adj=1)
matplot(-lp, cbind(qnU6, qnAsy) - qnorm(lp, lower.tail=TRUE, log.p=TRUE),
        pch=1:7, cex = .5, xaxt = "n", # and use eaxis() instead
        main = "absolute Error of qnorm() approximations", type = "b", log = "x")
sfsmisc::eaxis(1, sub10=2); p.ver()
legend("bottom", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:6, lty=1:5, pch=1:7, pt.cex=.5)

## If you look at the numbers, in versions of R <= 4.2.x,
## qnorm() is *worse* for large -lp than the higher order approximations
## ---> using qnormR() here:
absP <- function(re) pmax(abs(re), 2e-17) # not zero, so log-scale "shows" it
qnT <- qnormR(lp, lower.tail=TRUE, log.p=TRUE, version="2022") # ~= TRUE qnorm()
matplot(-lp, absP(cbind(qnU6, qnAsy) / qnT - 1),
        ylim = c(2e-17, .01), xaxt = "n", yaxt = "n", col=1:7, lty=1:7,
        main = "relative |Error| of qnorm() approximations", type = "l", log = "xy")
abline(h = .Machine$double.eps * c(1/2, 1, 2), col=adjustcolor("bisque",3/4),
       lty=c(5,2,5), lwd=c(1,3,1))
sfsmisc::eaxis(1, sub10 = 2, nintLog=20)
sfsmisc::eaxis(2, sub10 = c(-3, 2), nintLog=16)
mtext("qnT <- qnormR(*, version=\"2022\")", cex=0.9, adj=1)# ; p.ver()
legend("topright", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:7, lty=1:7, cex = 0.8)


###=== Optimal cut points / regions for different approximation orders k =================

## Zoom into each each cut-point region :
p.qnormAsy2 <- function(r0, k, # use k-1 and k in region around r0
                        n = 2048, verbose=TRUE, ylim = c(-1,1) * 2.5e-16,
                        rr = seq(r0 * 0.5, r0 * 1.25, length = n), ...)
{
  stopifnot(is.numeric(rr), !is.unsorted(rr), # the initial 'r'
            length(k) == 1L, is.numeric(k), k == as.integer(k), k >= 1)
  k.s <- (k-1L):k; nks <- paste0("k=", k.s)
  if(missing(r0)) r0 <- quantile(rr, 2/3)# allow specifying rr instead of r0
  if(verbose) cat("Around r0 =", r0,";  k =", deparse(k.s), "\n")
  lp <- (-rr^2) # = -r^2 = -s  <==> rr = sqrt(- lp)
  q. <- qnormR(lp, lower.tail=FALSE, log.p=TRUE, version="2022-08")# *not* depending on R ver!
  pq <- pnorm(q., lower.tail=FALSE, log.p=TRUE) # ~= lp
  ## the arg of pnorm() is the true qnorm(pq, ..) == q.  by construction
  ## cbind(rr, lp, q., pq)
  r <- sqrt(- pq)
  stopifnot(all.equal(rr, r, tol=1e-15))
  qnAsy <- sapply(setNames(k.s, nks), function(ord)
                  qnormAsymp(pq, lower.tail=FALSE, log.p=TRUE, order=ord))
  relE <- qnAsy / q. - 1
  m <- cbind(r, pq, relE)
  if(verbose) {
    print(head(m, 9)); for(j in 1:2) cat(" ..........\n")
    print(tail(m, 4))
  }
  ## matplot(r, relE, type = "b", main = paste("around r0 = ", r0))
  matplot(r, relE, type = "l", ylim = ylim,
     main = paste("Relative error of qnormAsymp(*, k) around r0 = ", r0,
                  "for  k =", deparse(k.s)),
     xlab = quote(r == sqrt(-log(p))), ...)
  legend("topleft", nks, col=1:2, lty=1:2, bty="n", lwd=2)
  for(j in seq_along(k.s))
    lines(smooth.spline(r, relE[,j]), col=adjustcolor(j, 2/3), lwd=4, lty=2)
  cc <- "blue2"; lab <- substitute(r[0] == R, list(R = r0))
  abline(v  = r0, lty=2, lwd=2, col=cc)
  axis(3, at= r0, labels=lab, col=cc, col.axis=cc, line=-1)
  abline(h = (-1:1)*.Machine$double.eps, lty=c(3,1,3),
         col=c("green3", "gray", "tan2"))
  invisible(cbind(r = r, qn = q., relE))
}

r0 <- c(27, 55, 109, 840, 36000, 6.4e8) # <--> in ../R/norm_f.R {and R's qnorm.c eventually}
## use k =   5   4    3    2      1       0    e.g.  k = 0  good for r >= 6.4e8
for(ir in 2:length(r0)) {
  p.qnormAsy2(r0[ir], k = 5 +2-ir) # k = 5, 4, ..
  if(interactive() && ir < length(r0)) {
       cat("[Enter] to continue: "); cat(readLines(stdin(), n=1), "\n") }
}

}
\keyword{distribution}
