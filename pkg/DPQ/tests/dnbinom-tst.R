#### Testing  1) dbinom_raw(), dnbinomR() and dnbinom.mu()
####          2) log1pmx(), logcf() etc
require(DPQ)
source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
## -> showProc.time(), assertError()

(doExtras <- DPQ:::doExtras() && !grepl("valgrind", R.home()))

if(!dev.interactive(orNone=TRUE)) pdf("wienergerm-accuracy.pdf")


### 1. Testing  dbinom_raw(), dnbinomR() and dnbinom.mu() >>> ../R/dbinom-nbinom.R <<<
### ----------   ../man/dbinom_raw.Rd & ../man/dnbinomR.Rd

## "FIXME:" use  sfsmisc :: relErrV() already here

###  dbinom() vs  dbinom.raw() :

for(n in 1:20) {
    cat("n=",n," ")
    for(x in 0:n)
        cat(".")
        for(p in c(0, .1, .5, .8, 1)) {
            stopifnot(all.equal(dbinom_raw(x, n, p, q=1-p, log=FALSE),
                                dbinom    (x, n, p,        log=FALSE)),
                      all.equal(dbinom_raw(x, n, p, q=1-p, log =TRUE),
                                dbinom    (x, n, p,        log =TRUE)))
    }
    cat("\n")
}
showProc.time()

###  dnbinom*() :
stopifnot(exprs = {
    dnbinomR(0, 1, 1) == 1
})

### exploring 'eps' == "true" tests must be done with  Rmpfr !!

### 2. Testing  log1pmx(), logcf() etc
### ----------

### 2a:  logcf()
##  ==   =======
x <- c((-20:3)/4, (25:31)/32) # close (but not too close) to upper bound 1

(lC   <- logcf (x, i=2, d=3, eps=1e-9))
 lCt  <- logcf (x, i=2, d=3, eps=1e-9, trace=TRUE) ; stopifnot(identical(lCt, lC))
(lR   <- logcfR(x, i=2, d=3, eps=1e-9))
all.equal(lC, lR, tol = 0) # to see if ..
stopifnot(all.equal(lC, lR, tol = 4e-16))
lRt  <- logcfR(x, i=2, d=3, eps=1e-9, trace=TRUE) ; stopifnot(identical(lRt, lR))
lRt2 <- logcfR(x, i=2, d=3, eps=1e-9, trace= 2)   ; stopifnot(identical(lRt2,lR))

lR.  <- logcfR.(x, i=2, d=3, eps=1e-9)
lR.t <- logcfR.(x, i=2, d=3, eps=1e-9, trace=TRUE) ; stopifnot(identical(lR.t, lR.))

all.equal(lC, lR., tol = 0) # TRUE !! (every where ?)
all.equal(lR, lR., tol = 0) # TRUE !!   "    "
stopifnot(all.equal(lC, lR., tol = 1e-15))
## (even though they used eps=1e-9 .. i.e., are not *so* accurate)
showProc.time()

##--- now with improved logcfR.() {<< will become the new logcfR() at least for MPFR !}:

##require(Rmpfr) may be not, see if NS loading (via "::") is sufficient:
requireNamespace("Rmpfr") || quit("no")
##                -----      ----------
xM <- Rmpfr::mpfr(x, 512)
(ct.14 <- system.time(lR.14 <- logcfR.(xM, i=2, d=3, eps=1e-20, trace=TRUE))) # 0.55 sec
(ct14  <- system.time(lR14  <- logcfR (xM, i=2, d=3, eps=1e-20, trace=TRUE))) # 4    sec

all.equal(lR.14, lR14, tol=0) #  TRUE
identical(lR.14, lR14)        #  TRUE !! (not sure if on all platforms!)

SS <- function(ch, digits=7)
    sub(paste0("([0-9]{1,",digits,"})[0-9]*e"), "\\1e", ch)
## double prec <--> MPFR:    vvvv (same eps)
lR.9 <- logcfR.(xM, 2,3, eps=1e-9)
## show:
SS(Rmpfr::all.equal(Rmpfr::roundMpfr(lR.9, 64), lR, tol=0))# .. 5.1138e-16
stopifnot(Rmpfr::all.equal(lR.9,  lR  , tol=1e-15))
SS(Rmpfr::all.equal(lR.14, lR.9, tol=0))# .. 3.701..e-10
stopifnot(Rmpfr::all.equal(lR.14, lR.9, tol=1e-9))
showProc.time()


### 2b:  log1pmx()
##  ==   =========

## The first part is from original ../man/log1pmx.Rd

e <- if(doExtras) 2^-12 else 2^-8; by.p <- 1/(if(doExtras) 256 else 64)
xd <- c(seq(-1+e, 0+100*e, by=e), seq(by.p, 5, by=by.p)) # length 676 or 5476 if do.X.
plot(xd, log1pmx(xd), type="l", col=2, main = "log1pmx(x)")
abline(h=0, v=-1:0, lty=3)

xM <- Rmpfr::mpfr(xd, 512)
asNumeric <- Rmpfr::asNumeric
## for MPFR numbers, really need more than tol_logcf = eps = 1e-14 (default)
if(doExtras) {
  print( system.time( # the logcfR() completely unvectorized is *slow*
        lg1pM <- log1pmx(xM, tol_logcf = 1e-25, eps2 = 1e-4, trace=TRUE)
  )) # 2 sec [new logcfR.()] if(doExtras) 0.43s otherwise, but ??s (!) on winbuilder!
  lines(xd, asNumeric(lg1pM), col=adjustcolor(4, 3/4), lwd=2)
}
lg1pdR <- log1pmx(xd, trace=TRUE, logCF=logcfR)
lines(xd, lg1pdR, col=adjustcolor(5, 3/4), lwd=3)
##
lg1pM. <- log1p(xM) - xM
xM2k <- Rmpfr::mpfr(xd, 2048) # even more accurate
lg1pM2k <- log1p(xM2k) - xM2k
relErrV <- sfsmisc::relErrV
reE00 <- asNumeric( relErrV(lg1pM2k, lg1pM.) )
print(signif(range(reE00)), 2) # [-1.5e-151, 4.8e-151] -- 512 bits is perfect
if(doExtras) { ## the error of the log1pmx() "algorithm" even for "perfect" mpfr-accuracy:
  rE.log1pm <- asNumeric(relErrV(lg1pM2k, lg1pM))
  print(signif(range(  rE.log1pm ), 2)) # -1.1e-26 7.9e-28
}
re <- asNumeric(relErrV(lg1pM., log1pmx(xd)))
showProc.time()

## MM: From around here, "move" to vignette  <<<<<<<<<<<<<<<<<<<<<<<<< ../vignettes/log1pmx-etc.Rnw <<<<<<

plot(xd, re, type="b", cex=1/4,
     xlab=quote(x), main="relErrV(log1p(xM) - xM, log1pmx(x)),  xM <- mpfr(x, 2048)")
abline(h = (-2:2)*2^-52, lty=2, col=adjustcolor("gray20", 1/2))
xl <- -0.84; xr <- -0.4 ##  the "relevant x-range":
colR <- adjustcolor("tomato", 1/2)
abline(v= c(xl,xr), col=colR, lty=2, lwd=2)
axis(1,at=c(xl,xr), labels=NULL, col.axis="tomato", col=colR, lwd=2, cex.axis=3/4)
text(xl, adj=c(0,-1/2), par("usr")[3], "zoom in", col=colR)
showProc.time()

# only "relevant"  x-range: values inside (xl, xr)
iN <- xl < xd & xd < xr
plot(xd[iN], re[iN], type="b", cex=1/4, main="rel.Error of log1pmx(x)",
     xlab=quote(x), ylab=quote(rE(x)), sub = "(in relevant x-range)")
abline(h = (-2:2)*2^-52, lty=2, col=adjustcolor("gray30", 1/2))
abline(v = c(xl,xr), col=colR, lty=2, lwd=2)

plot(xd[iN], abs(re[iN]), type="b", cex=1/2, log="y",
     main = "| relErr( log1pmx(x) ) |  {via 'Rmpfr'}",
     ylab=quote(abs(rE(x))), xlab=quote(x),
     ylim = c(4e-17, max(abs(asNumeric(re)[iN]))))
abline(h = c(1:2,4)*2^-52, lty=2, col=adjustcolor("gray20", 1/2))
mL1 <- eval(formals(log1pmx)$minL1)
abline(v = mL1, lwd=3, col=adjustcolor(2, 1/2))
axis(3, at=mL1, "minL1", col.axis=2, col=2)

mL1 <- -0.7; re2 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = mL1)))
lines(xd[iN], abs(re2), col=adjustcolor(4, 1/2), lwd=2)
abline(v = mL1, lwd=3,  col=adjustcolor(4, 2/3), lty=3)
axis(3, line=-1, at=mL1, sprintf("mL1 = %g",mL1), col.axis=4, col=4)
R <- xd[iN] >= mL1; stopifnot(re2[R] == re[iN][R])

mL1 <- -0.66; re3   <- asNumeric(relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = mL1)))
              re3.2 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = mL1, eps2 = .02)))
       print(max(abs(re3))) == # 4.818203e-16
             max(abs(re3.2)) #  TRUE !  of course, eps2 does not matter as long as it is far from our x-range !!
xd[iN][which.max(abs(re3))] #  -0.5710449
lines(xd[iN], abs(re3), col=adjustcolor(6, 1/3), lwd=2)
abline(v = mL1, lwd=3,  col=adjustcolor(6, 2/3), lty=3)
axis(3, line=-2, at=mL1, sprintf("mL1 = %g",mL1), col.axis=6, col=6)
lines(lowess(xd[iN], abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
lines(lowess(xd[iN], abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
lines(lowess(xd[iN], abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
R <- xd[iN] > mL1;  stopifnot(re3[R] == re2[R])
## MM: for vignette, stop here?
## ---> since  max(|re4|) {below} == max(|re3|)

mL1 <- -0.64; re4 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = mL1)))
max(abs(re4)) == max(abs(re3)) # TRUE !
lines(xd[iN], abs(re4), col=adjustcolor(7, 1/3), lwd=2)
abline(v = mL1, lwd=3,  col=adjustcolor(7, 2/3), lty=3)
axis(3, line=-2, at=mL1, sprintf("mL1 = %g",mL1), col.axis=7, col=7)
lines(lowess(xd[iN], abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
lines(lowess(xd[iN], abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
lines(lowess(xd[iN], abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
lines(lowess(xd[iN], abs(re4),    f=1/50), col=adjustcolor(7, 1/2), lwd=6)
R <- xd[iN] > mL1;  stopifnot(re4[R] == re3[R])

re5 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(xd[iN], minL1 = -0.6)))
max(abs(re5)) == max(abs(re3)) # TRUE !.. still not better
lines(xd[iN], abs(re5), col=adjustcolor(8, 1/3), lwd=2)
abline(v = -0.6, lwd=3, col=adjustcolor(8, 2/3), lty=3)
axis(3, line=-1, at=-0.6, sprintf("mL1 = %g",mL1), col.axis=8, col=8)
lines(lowess(xd[iN], abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
lines(lowess(xd[iN], abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
lines(lowess(xd[iN], abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
lines(lowess(xd[iN], abs(re4),    f=1/50), col=adjustcolor(7, 1/2), lwd=6)
lines(lowess(xd[iN], abs(re5),    f=1/50), col=adjustcolor(8, 1/2), lwd=6)
R <- xd[iN] > -0.6;  stopifnot(re5[R] == re4[R])
## -0.6  now is clearly too large, -0.7 was better  {MM: why the
showProc.time()


##--- separate: is  eps2 = 0.01  "optimal" / "good enough" ??

k. <- if(doExtras) 3 else 0
str(xS <- unique(sort(c(seq(-1/32, 1/32, by=2^-(12+k.)),
                        seq(-1,1, by=2^-(7+k.))/2^(4+k.))))) ## 257 or 3841 values
plot(xS, log1pmx(xS), type="l", col=2, main = "log1pmx(x)")
abline(h=0, v = c(-.01, .01), lwd=3, col=adjustcolor(7, 2/3), lty=3)
## nothing visible at +- .01

xSm <- Rmpfr::mpfr(xS, 1024)
xSM <- Rmpfr::mpfr(xS, 8192)
lg1pM. <- log1p(xSm) - xSm ## the direct formally using many bits ==> not much cancellation
summary(abs(asNumeric(relErrV(log1p(xSM) - xSM, lg1pM.))))# <= 4.1e-304 : 1024 is enough!
reS <- asNumeric(relErrV(lg1pM., log1pmx(xS)))
summary(reS) ## Lnx 64b gcc 11.2.1, 20210728:
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
## -2.625e-16 -5.246e-17  1.973e-18  1.211e-18  5.286e-17  2.688e-16
plot(xS, abs(reS), type="b", cex=1/2, log="y",
     main = "| relErr( log1pmx(<small x>) ) |  {via 'Rmpfr'}"
   , ylim = c(4e-17, max(abs(asNumeric(reS))))
   , panel.last= abline(h=0, v = c(-.01, .01), lwd=3, col=adjustcolor(7, 2/3), lty=3))
## no problem visible at  +/- eps2  =  +/- 0.01
stopifnot(abs(reS) < 1e-15, mean(abs(reS)) < 6e-16) #  seen 2.688e-16, 6.35e-17
showProc.time()

## in spite of the above; ... what's the approximation error of accurate log1pmx() ?
p.Err <- function(eps2 = .01, tol_logcf = 1e-17, xM=xSm) {
    stopifnot(eps2 > 0, tol_logcf > 0,
              inherits(xM, "mpfr"), (prec <- Rmpfr::getPrec(xM)) >= 64)
    lg1pM <- log1p(xM) - xM # direct formula suffering from cancellation
    lg1pM.<- log1pmx(xM, tol_logcf=tol_logcf, eps2=eps2)
    sys.call()
    reM <- asNumeric(relErrV(lg1pM., lg1pM))
    plot(asNumeric(xM), abs(reM), type="b", cex=1/2, log="y", xlab = quote(x),
         ylim = quantile(abs(reM), c(1/20,1), names=FALSE),
         sub = paste("min(getPrec(x)) =", min(prec)),
         main = sprintf("relative error of R log1pmx(eps2=%g, tol_logcf=%g)", eps2, tol_logcf))
    abline(h=0, v = eps2*c(-1,1), lwd=3, col=adjustcolor(7, 2/3), lty=3)
    invisible(reM)
}

### --- for tol_logcf = 1e-17 ---------------------------------------------
reSM    <- p.Err()     #  aha! Here,  eps2 = .01  was definitely too large!
reSM001 <- p.Err(.001) # still too large
xM2 <- Rmpfr::mpfr(seq(-1,1, length=1+1024)/2^8, 1024)
reSM3_4 <- p.Err(3e-4,    xM=xM2)# still too large
reSM294 <- p.Err(2.9e-4,  xM=xM2)# still
reSM2854<- p.Err(2.85e-4, xM=xM2)# still (and practically identical to 2.9e-4)
## cbind(reSM294, reSM2854)[reSM294 != reSM2854 , ]
##            reSM294      reSM2854
## [1,]  2.572818e-36  2.988004e-46
## [2,] -2.565617e-36 -2.977915e-46
reSM284 <- p.Err(2.8e-4,  xM=xM2)# now too small -- bad!
reSM274 <- p.Err(2.7e-4,  xM=xM2)#   (ditto) too small
showProc.time()

### --- for tol_logcf = 1e-14 == R's default!
reSM     <- p.Err(       tol = 1e-14) ## eps2=.01 is *still* too large
reSM002  <- p.Err(.002,  tol = 1e-14) # still too large
reSM0017 <- p.Err(.0017, tol = 1e-14) # still ..
reSM00165<- p.Err(.00165,tol = 1e-14) # still..
reSM00163<- p.Err(.00163,tol = 1e-14) # still..  <<<<<<<<<<<<< SHOULD WE CHANGE R TO USE THIS?
reSM00162<- p.Err(.00162,tol = 1e-14) # already (very barely!) too small
reSM0016 <- p.Err(.0016, tol = 1e-14) # already (barely!) too small
showProc.time()

summary(warnings())
