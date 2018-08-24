### Orig. = R script /u/maechler/R/MM/NUMERICS/dpq-functions/wienergerm-accuracy.R

library(DPQ)

if(!dev.interactive(orNone=TRUE)) pdf("wienergerm-accuracy.pdf")

nuS <- c(1:3, 10, 20, 30,50, 100, 200, 1000, 5000, 1e4)
ncpS <- c(0, .1, 1, 3, 10, 30, 100, 300, 1000, 10000, 1e5)
facs <- c(outer(c(1,2,5), 10^c(-1:3)))
r <- array(NA, dim = c(length(facs), 3, length(nuS), length(ncpS)),
           dimnames= list(x=NULL, c("p.","pW","pP"),
                          df = formatC(nuS), ncp= formatC(ncpS)))
for(df in nuS) {
    for(ncp in ncpS) {
        m <- df + ncp
        s <- sqrt(2*(df + 2*ncp))
        x <- m + s * facs
        r[,, formatC(df), formatC(ncp)] <-
            cbind(pchisq     (x, df=df, ncp=ncp),
                  pchisqW    (x, df=df, ncp=ncp),
                  pnchisq.Pea(x, df=df, ncp=ncp))
    }
    cat("done df=",formatC(df),"\n")
}

p.pchisq <- function(x, df, ncp, log = "x",
                     relErr = FALSE, diff = FALSE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 18 Mar 2004, 17:04

    Px <- pchisq(x, df=df, ncp=ncp)
    pW  <- pchisqW(x, df=df, ncp=ncp) # "Wiener2"
    pW1 <- pchisqW(x, df=df, ncp=ncp, variant = 'f')
    pP1 <- pnchisq.Pea(x, df=df, ncp=ncp)
    pP2 <- pnchisq.Pat(x, df=df, ncp=ncp)
    Pmat <- cbind(pW,pW1,pP1,pP2)

    tit <- paste("pchisq[Appr](x, df=",formatC(df),", ncp=",formatC(ncp),")")
    if(diff)
        tit <- paste(tit," -  pchisq[R](.)")
    else if( relErr)
        tit <- paste("rel.Error of", tit," to pchisq[R](.)")
    log.y <- length(grep("y", log))
    col <- 2:5
    lty <- 1:4
    noR <- diff || relErr
    if(noR) {
         Pmat <- Pmat - Px
         if(relErr)
             Pmat <- abs(Pmat / Px)
         else if(log.y) {
             Pmat <- abs(Pmat)
             if(log.y) tit <- paste("|", tit, "|")
         }
    }
    else {
        Pmat <- cbind(Px, Pmat)
        col <- c(1,col)
        lty <- c(1,lty)
    }
    if(log.y) op <- options(warn=-1)# no "log warning"
    matplot (x, Pmat, log=log, type='l', main = tit,
             ylab = '', ## ylim = rrange(Pmat, range = 0.1),
             lty = lty, col = col)
    if(log.y) op <- options(op)
    u <- par("usr")
    log.x <- length(grep("x", log))
    if(log.x) u[1:2] <- 10^u[1:2]
    if(log.y) u[3:4] <- 10^u[3:4]
    legend(u[1],u[4],
           c(if(!noR)"R", "Wiener2", "Wiener1", "Pearson","Patnaik"),
           lty = lty, col = col)
}

### all these for (0.1, 0.1):
x <- 2^seq(-10,3, length=501)
p.pchisq(x, df=0.1, ncp=0.1)
p.pchisq(2^seq(-12,6, length=2001), df=0.1, ncp=0.1)
p.pchisq(2^seq(-1, 6, length=2001), df=0.1, ncp=0.1)
p.pchisq(2^seq(-1, 6, length=2001), df=0.1, ncp=0.1, diff = TRUE)

p.pchisq(2^seq(-1, 10, length=2001), df=0.1, ncp=0.1, diff = TRUE)
p.pchisq(2^seq(-1, 10, length=2001), df=0.1, ncp=0.1, diff = TRUE,log="xy")
p.pchisq(2^seq(-5,  8, length=2001), df=0.1, ncp=0.1, diff = TRUE,log="xy")
p.pchisq(2^seq(-5,  8, length=2001), df=0.1, ncp=0.1, relE = TRUE,log="xy")

## (1, 1)
p.pchisq(2^seq(-14,8, length=2001), df=1, ncp= 1)
p.pchisq(2^seq(-14,8, length=2001), df=1, ncp= 1,diff=TRUE)
p.pchisq(2^seq(-20,9, length=2001), df=1, ncp= 1,diff=TRUE,log="xy")
p.pchisq(2^seq(-20,9, length=2001), df=1, ncp= 1,relE=TRUE,log="xy")

## (1, 10)
p.pchisq(2^seq(-14,8, length=2001), df=1, ncp= 10,diff=TRUE,log="xy")
## (10, 1)
p.pchisq(2^seq(-14,8, length=2001), df=10, ncp= 1,diff=TRUE,log="xy")

## (10, 10)
p.pchisq(2^seq(-14,8, length=2001), df=10, ncp= 10,diff=TRUE,log="xy")
p.pchisq(2^seq(-1,10, length=2001), df=10, ncp= 10,diff=TRUE,log="xy")

## (10, 100)
p.pchisq(2^seq(-14,8, length=2001), df=10, ncp= 100,diff=TRUE,log="xy")
p.pchisq(2^seq(-1,10, length=2001), df=10, ncp= 100,diff=TRUE,log="xy")
## "good" news: Pearson is "nice" when Wiener fails

## (10, 1e4)
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, log='')
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, diff=TRUE)
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, diff=TRUE, log='y')
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, relE=TRUE, log='y')

## (1000, 1)
p.pchisq(seq(800, 1200, length=2001), df=1000, ncp= 1, log='')

p.pchisq(seq(800, 1200, length=2001), df=1000, ncp= 1, diff=TRUE, log='y')
p.pchisq(seq(800, 1200, length=2001), df=1000, ncp= 1, relE=TRUE, log='y')
## interesting: Pearson really better than Wiener in a whole region
## x in [900, 1100] ...
