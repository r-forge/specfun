### R code from vignette source '..../DPQ/vignettes/qnorm-asymp.Rnw'

###################################################
### code chunk number 26: p.qnormAsy2-def
###################################################
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
  pq <- pnorm (q., lower.tail=FALSE, log.p=TRUE) # ~= lp
  ## the arg of pnorm() is the true qnorm(pq, ..) == q.  by construction
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
  legend("topleft", nks, horiz = TRUE, col=1:2, lty=1:2, bty="n", lwd=2)
  for(j in seq_along(k.s))
    lines(smooth.spline(r, relE[,j]), col=adjustcolor(j, 2/3), lwd=4, lty="6132")
  cc <- "blue2"; lab <- substitute(r[0] == R, list(R = r0))
  abline(v  = r0, lty=2, lwd=2, col=cc)
  axis(3, at= r0, labels=lab, col=cc, col.axis=cc, line=-1)
  abline(h = (-1:1)*.Machine$double.eps, lty=c(3,1,3),
         col=c("green3", "gray", "tan2"))
  invisible(cbind(r = r, qn = q., relE))
}

r0 <- c(27, 55, 109, 840, 36000, 6.4e8) # <-- cutoffs  <--> in ../R/norm_f.R
# use k =  5   4    3    2      1       0    e.g.  k = 0  good for r >= 6.4e8
if(FALSE) # not here, but typically from caller
  for(ir in 2:length(r0))
    p.qnormAsy2(r0[ir], k = 5 +2-ir, verbose=FALSE, cex.main = .90)

