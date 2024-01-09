#### Testing  stirlerr(), bd0(), ebd0(), dpois_raw(), ...
#### ===============================================

require(DPQ)
for(pkg in c("Rmpfr", "DPQmpfr"))
    if(!requireNamespace(pkg)) {
        cat("no CRAN package", sQuote(pkg), " ---> no tests here.\n")
        q("no")
    }
require("Rmpfr")

source(system.file(package="DPQ", "test-tools.R", mustWork=TRUE))
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS()
##_ options(conflicts.policy = list(depends.ok=TRUE, error=FALSE, warn=FALSE))
require(sfsmisc) # masking  'list_' *and* gmp's factorize(), is.whole()
##_ options(conflicts.policy = NULL)

showProc.time()

cutoffs <- c(15,35,80,500) # cut points, n=*, in the stirlerr() "algorithm"
##
n <- c(seq(1,15, by=1/4),seq(16, 25, by=1/2), 26:30, seq(32,50, by=2), seq(55,1000, by=5),
       20*c(51:99), 50*(40:80), 150*(27:48), 500*(15:20))
st.n <- stirlerr(n)# rather use.halves=TRUE, just here , use.halves=FALSE)
plot(st.n ~ n, log="xy", type="b") ## looks good now (straight line descending {in log-log !}
nM <- mpfr(n, 2048)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose
all.equal(asNumeric(st.nM), st.n)# TRUE
all.equal(st.nM, as(st.n,"mpfr"))# .. difference: 3.381400e-14 was 1.05884.........e-15
all.equal(roundMpfr(st.nM, 64), as(st.n,"mpfr"), tolerance=1e-16)# (ditto)

## --- Look at the direct formula -- why is it not good for n ~= 5  ?
##
## Preliminary Conclusions :
## 1. it's not clear why
## 2. lgamma1p(n) does really not help much compared to lgamma(n+1) --- surprisingly !!

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param n
##' @return
##' @author Martin Maechler
relE.lgam1 <-  function(n) {
    nM <- mpfr(n, 1024)
    st  <- lgamma(n +1)  - (n +0.5)*log(n)  + n  - log(2*pi)/2
    st. <- lgamma1p(n)   - (n +0.5)*log(n)  + n  - log(2*pi)/2
    stM <- lgamma(nM+1)  - (nM+0.5)*log(nM) + nM - log(2*Const("pi", 256))/2
    ## stM <- roundMpfr(stM, 128)
    cbind(dir = asNumeric(relErrV(stM, st)),
          l1p = asNumeric(relErrV(stM, st.)))
}

n <- 2^-(1000:1)
relEx <- relE.lgam1(n)

cols <- c("gray30", adjustcolor(2, 1/2))
matplot(n, relEx, type = "l", log="x", col=cols, lty=1, lwd = c(1,3))
## very interesting: there are different intervals  <---> log(n) pattern !!
## -- but very small difference, only for n >~= 1/1000  but not before
abline(h = c(-4,-2:2, 4)*2^-53, lty=1, col="gray")

## Check  log() only :
plot(n, asNumeric(relErrV(log(mpfr(n, 256)), log(n))),
     log="x", type="l") ## ===> indeed --- log(n) approximation pattern !!

## really for the very small n, all is dominated by -(n+0.5)*log(n);  and lgamma1p() is unnecessary!
i <- 1:20; ni <- n[i]
lgamma1p(ni)
- (ni +0.5)*log(ni)  + ni

## much less extreme:
n2 <- lseq(2^-12, 1/2, length=1000)
relE2 <- relE.lgam1(n2)

cols <- c("gray30", adjustcolor(2, 1/2))
matplot(n2, relE2, type = "l", log="x", col=cols, lty=1, lwd = c(1,3))
abline(h = (-4:4)*2^-53, lty=1, col="gray")

matplot(n2, pmax(abs(relE2), 1e-17), type = "l", log="xy", col=cols, lty=1, lwd = c(1,3))
abline(h = (1:4)*2^-53, lty=1, col="gray")

## which is better?  ... "random difference"
d.absrelE <- abs(relE2[,"dir"]) - abs(relE2[,"l1p"])
plot (n2, d.absrelE, type = "l", log="x", # no clear picture ...
      main = "|relE_dir|  -  |relE_l1p|", axes=FALSE, frame.plot=TRUE)
eaxis(1, sub10=c(-2,1)); eaxis(2); axis(3, at=max(n2)); abline(v = max(n2), lty=3, col="gray")
## 'l1p' very slightly better:
lines(n, smooth.spline(d.absrelE, df=12)$y, lwd=3, col=2)

## not really small n at all == here see, how "bad" the direct formula gets for 1 < n < 10 or so
n3 <- lseq(2^-14, 2^2, length=800)
relE3 <- relE.lgam1(n3)

matplot(n3, relE3, type = "l", log="x", col=cols, lty=1, lwd = c(1,3))

matplot(n3, pmax(abs(relE3), 1e-17), type = "l", log="xy", col=cols, lty=1, lwd = c(1,3), xaxt="n")
abline(h = (1:4)*2^-53, lty=1, col="gray"); eaxis(1, sub10=c(-2,3))

## very small difference:
ll3.1 <- lowess(log(n3), log(abs(relE3[,1])), f= 0.1)
ll3.2 <- lowess(log(n3), log(abs(relE3[,2])), f= 0.1)

with(ll3.1, lines(exp(x), exp(y), col="darkblue", lwd=3))
with(ll3.2, lines(exp(x), exp(y), col="red3",     lwd=3))
## ==> ok:  lgamma1p(.) is very slightly better !!!




##' Very revealing plot showing the *relative* approximation error of stirlerr(<dblprec>)
##'
p.stirlerrDev <- function(n, precBits=2048,
                          stnM = stirlerr(mpfr(n, precBits), use.halves=use.halves, verbose=verbose),
                          abs = FALSE,
                          ## cut points, n=*, in the stirlerr() algorithm :
                          scheme = c("R3", "R4.x"),
                          cutoffs = switch(match.arg(scheme),
                                           R3 = c(15, 35, 80, 500),
                                           R4.x = c(7.5, 8.5, 10.625, 12.125, 20, 26, 55, 200, 3300)),
                          use.halves = missing(cutoffs),
                          verbose = getOption("verbose"),
                          type = "b", cex = 1,
                          col = adjustcolor(1, 3/4), colnB = adjustcolor("orange4", 1/3),
                          log = if(abs) "xy" else "x",
                          xlim=NULL, ylim = if(abs) c(8e-18, max(abs(N(relE)))))
{
    op <- par(las = 1, mgp=c(2, 0.6, 0))
    on.exit(par(op))
    st <- stirlerr(n, cutoffs=cutoffs, use.halves=use.halves)
    relE <- relErrV(stnM, st) # eps0 = .Machine$double.xmin
    N <- asNumeric
    form <- if(abs) abs(N(relE)) ~ n else N(relE) ~ n
    plot(form, log=log, type=type, cex=cex, col=col, xlim=xlim, ylim=ylim,
         ylab = quote(relErrV(stM, st)), axes=FALSE, frame.plot=TRUE,
         main = sprintf("stirlerr(n, cutoffs) rel.error [wrt stirlerr(Rmpfr::mpfr(n, %d))]",
                        precBits))
    eaxis(1, sub10=3)
    eaxis(2)
    mtext(paste("cutoffs =", deparse(cutoffs)))
    ylog <- par("ylog")
    if(ylog) {
        epsC <- c(1,2,4,8)*2^-52
        epsCxp <- expression(epsilon[C],2*epsilon[C], 4*epsilon[C], 8*epsilon[C])
    } else {
        epsC <- (-2:2)*2^-52
        epsCxp <- expression(-2*epsilon[C],-epsilon[C], 0, +epsilon[C], +2*epsilon[C])
    }
    dy <- diff(par("usr")[3:4])
    if(diff(range(if(ylog) log10(epsC) else epsC)) > dy/50) {
        lw <- rep(1/2, 5); lw[if(ylog) 1 else 3] <- 2
        abline(  h=epsC, lty=3, lwd=lw)
        axis(4, at=epsC, epsCxp, las=2, cex.axis = 3/4, mgp=c(3/4, 1/4, 0), tck=0)
    } else ## only x-axis
        abline(h=if(ylog) epsC else 0, lty=3, lwd=2)
    abline(v = cutoffs, col=colnB)
    axis(3, at=cutoffs, col=colnB, col.axis=colnB,
         labels = formatC(cutoffs, digits=3, width=1))
    invisible(relE)
} ## p.stirlerrDev()

do.pdf <- TRUE
do.pdf <- !dev.interactive(orNone = TRUE)
do.pdf
if(do.pdf)
    pdf("stirlerr-relErr_0.pdf", width=8, height=6)

showProc.time()

## store "expensive" stirlerr() result, and re-use many times below:
nM <- mpfr(n, 2048)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose

p.stirlerrDev(n=n, stnM=st.nM, use.halves = FALSE) # default cutoffs= c(15, 40, 85, 600)
## show the zoom-in region in next plot
yl2 <- 3e-14*c(-1,1)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

if(do.pdf) {
    dev.off() ; pdf("stirlerr-relErr_1.pdf", width=8, height=6)
}

## drop n < 5:
p.stirlerrDev(n=n, stnM=st.nM, xlim = c(7, max(n)), use.halves=FALSE) # default cutoffs= c(15, 40, 85, 600)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

## The first plot clearly shows we should do better:
## Current code is switching to less terms too early, loosing up to 2 decimals precision
p.stirlerrDev(n=n, stnM=st.nM, ylim = yl2, use.halves = FALSE)
p.stirlerrDev(n=n, stnM=st.nM, ylim = yl2, use.halves = TRUE)# exact at n/2 (n <= ..)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

showProc.time()


if(do.pdf) {
    dev.off(); pdf("stirlerr-relErr_6-fin-1.pdf")
}

### ~19.April 2021: "This is close to *the* solution" (but see 'cuts' below)
cuts <- c(7, 12, 20, 26, 60, 200, 3300)
##        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st.   <- stirlerr(n=n , cutoffs = cuts, verbose=TRUE)
st.nM <- stirlerr(n=nM, cutoffs = cuts, use.halves=FALSE) ## << on purpose
relE <- asNumeric(relErrV(st.nM, st.))
head(cbind(n, relE), 20)
## nice printout :
print(cbind(n       = format(n, drop0trailing = TRUE),
            stirlerr= format(st.,scientific=FALSE, digits=4),
            relErr  = signif(relE, 4))
      , quote=FALSE)

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts)

if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-2.pdf") }

## zoom in ==> {good for n >= 10}
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", ylim = 2e-15*c(-1,1),
              cutoffs = cuts)## old default cutoffs = c(15,35, 80, 500)

if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-3.pdf") }


##-- April 20: have more terms up to S10 in stirlerr() --> can use more cutoffs
n <- lseq(1/64, 5000, length=4096)
nM <- mpfr(n, 2048)
cuts <- c(5.4, 7.5, 8.5, 10.625, 12.125, 20, 26, 60, 200, 3300)
##        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (from below: 5.4 is too small..)
st.nM <- roundMpfr(stirlerr(nM, use.halves=FALSE, ## << on purpose;
                            verbose=TRUE), precBits = 128)
## NB: for x=xM <mpfr>; `cutoffs` are *not* used.
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts)

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, abs=TRUE)
## using exact values sferr_halves[] *instead* of MPFR ones
lines((0:30)/2, abs(stirlerr((0:30)/2, cutoffs=cuts, verbose=TRUE)/DPQ:::sferr_halves - 1), type="o", col=2,lwd=2)

if(do.pdf) { dev.off(); pdf("stirlerr-tst_others.pdf") }

## should we e.g., use interpolation spline through sfserr_halves[] for n <= 7.5
## -- doing the interpolation on the  log(1 - 12*x*stirlerr(x)) vs  log2(x)  scale -- maybe ?
stirL <- curve(1-12*x*stirlerr(x, cutoffs = cuts, verbose=TRUE), 1/64, 8, log="xy", n=2048)
## just need "true" values for x = 2^-(6,5,4,3,2) in addition to those we already have at x = 1/2, 1.5, 2, 2.5, ..., 7.5, 8
xM <- mpfr(stirL$x, 2048)
stilEM <- roundMpfr(1 - 12*xM*stirlerr(xM, verbose=TRUE), 128)
relEsml <- relErrV(stilEM, stirL$y)
plot(stirL$x, relEsml) # again: stirlerr() is "limited.." for ~ [2, 6]

## The function to interpolate:
plot(asNumeric(    (stilEM)) ~ x, data = stirL, type = "l", log="x")
plot(asNumeric(sqrt(stilEM)) ~ x, data = stirL, type = "l", log="x")
plot(asNumeric( log(stilEM)) ~ x, data = stirL, type = "l", log="x")

y <- asNumeric( log(stilEM))
x <- stirL$x
spl <- splinefun(log(x), y)

plot(asNumeric(log(stilEM)) ~ x, data = stirL, type = "l", log="x")
lines(x, spl(log(x)), col=2)
summary(rE <- relErrV(target = y, current = spl(log(x)))) # all 0 {of course: interpolation at *these* x}

ssmpl <- smooth.spline(log(x), y)
str(ssmpl$fit)# only 174 knots, 170 coefficients
methods(class="smooth.spline.fit") #-> ..
str(yp <- predict(ssmpl$fit, x=log(x)))
lines(x, yp$y, col=adjustcolor(3, 1/3), lwd=3) # looks fine
## but
summary(reSpl <- relErrV(target = y, current = yp$y))
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## -6.35e-10 -7.55e-11  1.10e-12  0.00e+00  7.32e-11  5.53e-10
## which is *NOT* good enough, of course ....

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim=c(-1,1)*4e-14)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim=c(-1,1)*1e-15)

st. <- stirlerr(n=n, cutoffs = cuts, verbose=TRUE)
relE <- asNumeric(relErrV(st.nM, st.))
plot(abs(relE) ~ n, type="o", pch=".", log="xy")

head(cbind(n, relE), 20)
## nice printout :
print(cbind(n       = format(n, drop0trailing = TRUE),
            stirlerr= format(st.,scientific=FALSE, digits=4),
            relErr  = signif(relE, 4))
      , quote=FALSE)

showProc.time()

##=== Really, dpois_raw()  and  dbinom_raw()  *both* use stirlerr(x)  for "all"  'x > 0'
##            ~~~~~~~~~~~       ~~~~~~~~~~~~             ===========              ===== !

## below, 7 "it's okay, but not perfect:" ===>  need more terms in stirlerr()  __or__ ??
## April 20: MM added more terms up to S10
x <- lseq(1/16, 7, length=2048)
system.time(stM <- DPQmpfr::stirlerrM(Rmpfr::mpfr(x,2048))) # 1.7 sec elapsed
plot(x, stirlerr(x, use.halves=FALSE) - stM, type="l", log="x", main="absolute Error")
plot(x,     stirlerr(x, use.halves=FALSE) / stM - 1,  type="l", log="x", main="relative Error")
plot(x, abs(stirlerr(x, use.halves=FALSE) / stM - 1), type="l", log="xy",main="|relative Error|")
abline(h=c(1,2,4)*.Machine$double.eps, lty=3)
## lgammacor() does *NOT* help, as it is  *designed*  for  x >= 10!
str(lgamc <- lgammacor(x)); table(is.nan(lgamc)) # *all* are NaN !
## maybe look at it for x >= 9 or so ?
##
## ==> Need another chebyshev() or rational-approx. for x in [.1, 7] or so !!

##=============> For now, see ../Misc/stirlerr-trms.R  <===============
##                            ~~~~~~~~~~~~~~~~~~~~~~~
showProc.time()
