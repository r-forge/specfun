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
           order, lower.tail = TRUE, log.p = FALSE)
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
  Martin Maechler (2022). Asymptotic Tail Formulas For Gaussian Quantiles;
  not yet finished.
}
\author{Martin Maechler}
\seealso{
  The upper tail approximations in Abramowitz & Stegun, in \pkg{DPQ}
  available as \code{qnormUappr()} and \code{\link{qnormUappr6}()},
  are less accurate than our \code{order >= 1} formulas in the tails.
}
\examples{%% MM: see also ./qnormUappr.Rd
lp <- -c(head(c(outer(c(5,2,1), 10^(18:1))), -2), 20:10, seq(9.75, 2, by = -1/8))
qnU6 <- qnormUappr6(lp=lp) # 'p' need not be specified if 'lp' is
qnAsy <- sapply(0:5, function(ord) qnormAsymp(lp=lp, lower.tail=FALSE, order=ord))
summary(warnings()) # 12 warnings ... NaNs produced
matplot(-lp, cbind(qnU6, qnAsy), type = "b", log = "x", pch=1:7)
legend("center", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:6, lty=1:6, pch=1:7)

matplot(-lp, cbind(qnU6, qnAsy) - qnorm(lp, lower.tail=TRUE, log.p=TRUE), pch=1:7,
        xaxt = "n", # and use eaxis() instead
        main = "absolute Error of qnorm() approximations", type = "b", log = "x")
sfsmisc::eaxis(1, sub10=2)
legend("bottom", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:6, lty=1:6, pch=1:7)

## If you look at the numbers, in versions of R <= 4.2.x,
## qnorm() is *worse* for large -lp than the higher order approximations

absP <- function(re) pmax(abs(re), 2e-17) # not zero, so log-scale "shows" it
matplot(-lp, absP(cbind(qnU6, qnAsy) / qnorm(lp, lower.tail=TRUE, log.p=TRUE) - 1),
        ylim = c(2e-17, .01), xaxt = "n", yaxt = "n", col=1:7, lty=1:7,
        main = "relative |Error| of qnorm() approximations", type = "l", log = "xy")
abline(h= .Machine$double.eps * c(1/2, 1, 2), col=adjustcolor("bisque",3/4), lty=c(5,2,5), lwd=c(1,3,1))
sfsmisc::eaxis(1, sub10 = 2, nintLog=20)
sfsmisc::eaxis(2, sub10 = c(-3, 2), nintLog=16)
legend("right", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:7, lty=1:7)

str(r <- c(1:63, 2^c(seq(6, 20, by=1/8), 21:60)))
## rounded / approximate cutoffs in r-scale for good qnorm() approx:
## original numbers %% from MM's NUMERICS/dpq-functions/qnorm-asymptotic.R
r0 <- c(27, 36, 61, 390, 14144)

## now used, after somewhat zooming into each cut-point region:
r0 <- c(57, 99, 840, 30000, 6e8)
for(ir in seq_along(r0)) {
  r. <- r0[ir]
  k <- 6 - ir # 5, 4, ..
  rr <- seq(r.*.80, r.*1.25, length = 2048)
  lp <- (-rr^2)
  q. <- qnorm(lp, lower.tail=FALSE, log.p=TRUE)
  pq <- pnorm(q., lower.tail=FALSE, log.p=TRUE) # ~= lp
  ##  <==>  true  qnorm(pq, ..) == q.  by construction
  ## cbind(rr, lp, q., pq)
  k.s <- (k-1):k; nks <- paste0("k=", k.s)
  cat("Around r0 =", r.,";  k =", deparse(k.s), "\n")
  qnAsy <- sapply(setNames(k.s, nks), function(ord)
                  qnormAsymp(pq, lower.tail=FALSE, log.p=TRUE, order=ord))
  relE <- qnAsy / q. - 1
  m <- cbind(rr, pq, relE)
  print(head(m, 9)); for(j in 1:2) cat(" ..........\n"); print(tail(m,4))
  ## matplot(rr, relE, type = "b", main = paste("around r0 = ", r.))
  matplot(rr, relE, type = "l",
     main = paste("Relative error qnormAsymp(*, k) around r0 = ", r., ";  k =", deparse(k.s)),
     xlab = quote(r == sqrt(-log(p))), ylim = c(-1,1)*2.5e-16)
  legend("topleft", nks, col=1:2, bty="n", lwd=2)
  for(j in seq_along(k.s))
    lines(smooth.spline(rr, relE[,j]), col=adjustcolor(j, 2/3), lwd=4, lty=2)
  abline(v = r., lty=2, col="gold")
  abline(h = (-1:1)*.Machine$double.eps, lty=c(3,1,3), col=c("green3", "gray", "tan2"))
  if(interactive() && r. != r0[length(r0)]) {
       cat("[Enter] to continue: ");  cat(readLines(stdin(), n=1), "\n") }
}

}
\keyword{distribution}
