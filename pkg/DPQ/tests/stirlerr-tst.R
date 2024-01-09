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
relErrV <- sfsmisc::relErrV

pks <- c("sfsmisc", "DPQ", "Rmpfr", "DPQmpfr")
sapply(lapply(setNames(,pks), packageVersion), format)

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

plot (n, abs(relE) - abs(relE.), type = "l", log="x") # no clear picture
lines(n, smooth.spline(abs(relE) - abs(relE.))$y, lwd=3, col=2)

## not really small n at all == here see, how "bad" the direct formula gets for 1 < n < 10 or so
n3 <- lseq(2^-14, 2^2, length=800)
relE3 <- relE.lgam1(n3)

matplot(n3, relE3, type = "l", log="x", col=cols, lty=1, lwd = c(1,3))

matplot(n3, pmax(abs(relE3), 1e-17), type = "l", log="xy", col=cols, lty=1, lwd = c(1,3))
abline(h = (1:4)*2^-53, lty=1, col="gray")

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
    sfsmisc::eaxis(1, sub10=3)
    sfsmisc::eaxis(2)
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
}

do.pdf <- TRUE
do.pdf <- !dev.interactive(orNone = TRUE)
do.pdf
if(do.pdf)
    pdf("stirlerr-relErr_0.pdf", width=8, height=6)

showProc.time()

p.stirlerrDev(n=n, stnM=st.nM, use.halves = FALSE) # default cutoffs= c(15, 40, 85, 600)
## show the zoom-in region in next plot
yl2 <- 3e-14*c(-1,1)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

if(do.pdf) {
    dev.off() ; pdf("stirlerr-relErr_1.pdf", width=8, height=6)
}

## drop n < 5:
p.stirlerrDev(n=n, stnM=st.nM, xlim = c(7, max(n)), use.halve=FALSE) # default cutoffs= c(15, 40, 85, 600)
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
st. <- stirlerr(n=n, cutoffs = cuts, verbose=TRUE)
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
n <- sfsmisc::lseq(1/64, 5000, length=4096)
nM <- mpfr(n, 2048)
st.nM <- roundMpfr(stirlerr(nM, use.halves=FALSE, ## << on purpose
                            verbose=TRUE), precBits = 128)

cuts <- c(5.4, 7.5, 8.5, 10.625, 12.125, 20, 26, 60, 200, 3300)
##        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (from below: 5.4 is too small..)
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
x <- sfsmisc::lseq(1/16, 7, length=2048)
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




###--------------- bd0()  & ebd0() ------------------------------------------------------


## ebd0 constants: the column sums of "bd0_scale":  log(n / 1024) for all these n
## ---- according to the *comments* in the C code -- so here I test that at least the *sums* are correct
bd0.n <- c(2048,2032,2016,2001,1986,1971,1956,1942,1928,1913,1900,1886,1872,1859,
           1846,1833,1820,1808,1796,1783,1771,1759,1748,1736,1725,1713,1702,1691,
           1680,1670,1659,1649,1638,1628,1618,1608,1598,1589,1579,1570,1560,1551,
           1542,1533,1524,1515,1507,1498,1489,1481,1473,1464,1456,1448,1440,1432,
           1425,1417,1409,1402,1394,1387,1380,1372,1365,1358,1351,1344,1337,1331,
           1324,1317,1311,1304,1298,1291,1285,1279,1273,1266,1260,1254,1248,1242,
           1237,1231,1225,1219,1214,1208,1202,1197,1192,1186,1181,1176,1170,1165,
           1160,1155,1150,1145,1140,1135,1130,1125,1120,1116,1111,1106,1101,1097,
           1092,1088,1083,1079,1074,1070,1066,1061,1057,1053,1049,1044,1040,1036,
           1032,1028,1024)

stopifnot(
    all.equal(bd0.n,
              1024 * exp(colSums(DPQ:::logf_mat)))
) ## on lynne (64-bit, Fedora 32, 2021) they are even *identical*
identical(bd0.n, 1024 * exp(colSums(DPQ:::logf_mat))) # amazingly to me

## Also, the numbers themselves decrease monotonely,
## their differences are close to, but *not* monotone:
diff(bd0.n) # -16 -16 -15 -15 -15 -15 -14 -14 -15 -13 -14 ...
            #                             ^^^^^^^^^^^^^^  (etc)

if(do.pdf) { dev.off(); pdf("diff-bd0_tab.pdf") }

plot(diff(bd0.n), type="b")
c2 <- adjustcolor(2, 1/2)
par(new=TRUE)
plot(diff(bd0.n, differences = 2), type="b", col=c2, axes=FALSE, ann=FALSE)
axis(4, at=-1:2, col=c2, col.axis=c2)

showProc.time()


## use functionality originally in ~/R/MM/NUMERICS/dpq-functions/15628-dpois_raw_accuracy.R
## now
require(Rmpfr)


## transition till DPQmpfr exports this *and* is on CRAN, to ease maintainer("DPQ")
##                     vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
if(file.exists(ff <- "~/R/Pkgs/DPQmpfr/R/dpoisEr.R")) withAutoprint({ #-------------

source(ff)
##-- ----- *or* move to vignette ../vignettes/log1pmx-etc.Rnw  <<<<<<<<<<<<<<<

##-------- small lambda --- is  dpois_simpl0() good ?

range(dpE40 <- dpoisEr(40.25, x=0:200))
## -4.401959e-16  3.645529e-16

## dpois_simpl0() uses "old" direct formula on original scale: factorial(x)
stopifnot(factorial(170) < Inf,
          factorial(171) == Inf)
xS <- 0:170 # the full range of "sensible" x values for dpois_simpl0
range(dpE40Sim <- dpoisEr(40.25, x=xS, dpoisFUN = dpois_simpl0))
## -1.299950e-13  1.118291e-13
p.dpoisEr(dpE40Sim)
## --> very suprising: errors are *very* small  up to  x <= 49, then in in the order of 1e-13
p.dpoisEr(dpE40Sim, ylim = c(-1,1)*4e-16)

## zoom in --- integer x only:
range(dpE40Sim2 <- dpoisEr(40.25, x=0:49, dpoisFUN = dpois_simpl0))
## -2.889661e-16  2.076597e-16
p.dpoisEr(dpE40Sim2) #   --- almost all in [-eps, +eps]

## zoom in and use non-integer x:
range(dpE40Sim2d <- dpoisEr(40.25, x=seq(0, 49, by=1/8), dpoisFUN = dpois_simpl0))
## [1] -3.861877e-14  3.080750e-14  == Oops !  blown up to
p.dpoisEr(dpE40Sim2d)

})

### MM: FIXME -- moved all this to an  Rmpfr  vignette:
###
" ~/R/D/R-forge/Rmpfr/pkg/vignettes/gamma-inaccuracy.Rnw "
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## and the R code actually to the R script:
" ~/R/D/R-forge/Rmpfr/pkg/vignettes/gamma-inaccuracy_src/plot-factErr-def.R "

##=======>  gamma(x) itself suffers from the fact that  exp(y) has a *large* relative error,
##          -------- when  |y| ~ 100 or so, more specifically, the
## relative error of   exp(y) =  |y| * {rel.err(y)} , since
##   exp(((1+ eps)*y) = exp(y) * exp(eps*y) >=  exp(y) (1 + eps*y)  and indeed,
## the inaccuracy of y (i.e. eps)  is blown up by a factor |y|  which is not small here!


## close to over-/underflow -------

### Large lambda == np == M -------

if(do.pdf) { dev.off(); pdf("stirlerr-bd0-ebd0.pdf") }

##-- TODO ----- *or* move to vignette ---> ../vi???????????

LL <- 1e20
dput(x1 <- 1e20 - 2e11) # 9.99999998e+19

(P1 <-         dpois     (x1,       LL)) # was 3.989455e-11; now 5.520993e-98
(P1m <- Rmpfr::dpois(mpfr(x1, 128), LL)) # 5.52099285934214335003128935..e-98
## However -- the ebd0() version
(P1e <- dpois_raw(x1, LL, version="ebd0_v1"))
## was 3.989455e-11, but now good!
stopifnot(exprs = {
    all.equal(P1 , 5.520992859342e-98, tol=1e-12)
    all.equal(P1e, P1, tol=1e-12)
    all.equal(P1m, P1, tol=1e-12)
})

options(digits=9)

## indeed:  regular  bd0()  works "ok"  --- but  ebd0() does non-sense !!
(bd.1 <- bd0(x1, LL, verbose=2))
## bd0(1e+20, 1e+20): T.series w/ 2 terms -> bd0=200
## [1] 200
(bd.1M <- bd0(x1, mpfr(LL, 128), verbose=2))
## bd0(1e+20, 1e+20): T.series w/ 3 terms -> bd0=200
## ---> 199.9999919413334091607468236761591740489
    asNumeric(bd.1 / bd.1M - 1)# -1.82e-17 -- suggests bd0() is really accurate here
stopifnot(abs(bd.1 / bd.1M - 1) < 3e-16,
          all.equal(199.999991941333, bd.1, tolerance=1e-14))

ebd0(x1, LL, verbose=TRUE)# fixed since  June 6, 2021

showProc.time()

### Large x  -- small  np == M ------------------------------------


mpfrPrec <- 1024
mpfrPrec <- 256

yy <-   bd0 (1e307, 10^(-5:1), verbose=TRUE)
yhl <- ebd0 (1e307, 10^(-5:1), verbose=TRUE)
yhlC<- ebd0C(1e307, 10^(-5:1)) ## newly returns data.frame
stopifnot(yy == Inf, yhl$yh == Inf, yhlC == yhl)
yM <- bd0(mpfr(1e307, mpfrPrec), 10^(-5:1))
roundMpfr(range(yM), 12) ##  7.0363e+309 7.1739e+309 -- *are* larger than  DBL_MAX


### Now  *BOTH*  x and lambda are large : ---------------------------------------
## (FIXME?? Small loss for ebd0, see below) <<< ???
##  is bd0(<mpfr>, *) really accurate ???
##  it uses it's own convergent series approximation for |x-np| < .. ????

x. <- 1e307
ebd0(x., 10^(300:308))

stopifnot(is.finite(Llam <- 2^(990:1024 - 1e-12)))

bd0ver <- function(x, np, mpfrPrec, chkVerb=TRUE, keepMpfr=FALSE) {
    stopifnot(length(mpfrPrec <- as.integer(mpfrPrec)) == 1,
              !is.na(mpfrPrec), mpfrPrec >= 64,
              x >= 0, np >= 0)
    yy   <-  bd0 (x, np)
    yhl  <- ebd0 (x, np)
    yhlC <- ebd0C(x, np)
    if(chkVerb) {
        yhl.  <- ebd0 (x, np, verbose=TRUE)
        yhlC. <- ebd0C(x, np, verbose=TRUE)
        stopifnot(identical(yhl., yhl),
                  identical(yhlC., yhlC))
    }
    epsC <- .Machine$double.eps
    aeq0 <- all.equal(yhl, yhlC, tol = 0)
    aeq4 <- all.equal(yhl, yhlC, tol = 4*epsC)
    if(!isTRUE(aeq4)) warning("the C and R versions of ebd0() differ:", aeq4)
    stopifnot(is.whole(yhl [["yh"]]),
              is.whole(yhlC[["yh"]]))
    yM  <-  bd0(mpfr(x, mpfrPrec),
                mpfr(np,mpfrPrec), verbose=chkVerb)# more accurate ! (?? always ??)
    relE <- relErrV(target = yM, # the mpfr one
                    cbind(ebd0 = yhl [["yh"]] + yhl [["yl"]],
                          ebd0C= yhlC[["yh"]] + yhlC[["yl"]],
                          bd0 = yy))
    relE <- structure(asNumeric(relE), dim=dim(relE), dimnames=dimnames(relE))
    ## return:
    list(x=x, np=np, bd0=yy, ebd0=yhl, ebd0C=yhlC,
         bd0M=if(keepMpfr) yM, # <- expensive
         aeq0=aeq0, aeq4=aeq4, relE = relE)
}

bd0v.8  <- bd0ver(x., Llam, mpfrPrec = 256)
bd0v.10 <- bd0ver(x., Llam, mpfrPrec = 1024)
stopifnot( all.equal(bd0v.8, bd0v.10, tol=0),
          bd0v.8$aeq0, # even tol=0 equality !
          bd0v.8$aeq4 )
## ==> 256 bit gives the *same* (asNumeric() - double-prec accuracy) as 1024 bits !
rm(bd0v.10)
showProc.time()

p.relE <- function(bd0v, dFac = if(max(np) >= 8e307) 1e10 else 1,
                   log = "x", type="b") {
    stopifnot(length(x <- bd0v$x) == 1 # for now
            , is.numeric(x), is.numeric(np <- bd0v$np), length(np) > 1
            , is.numeric(dFac), dFac > 0, length(dFac) == 1
            , is.matrix(relE <- bd0v$relE)
            , (k <- ncol(relE)) >= 1
            , sum(iOk <- local({ y <- bd0v$bd0; is.finite(y) & y != 0 })) > 1
              )
    ## */dFac : otherwise triggering axis() error
    ## log - axis(), 'at' creation, _LARGE_ range: invalid {xy}axp or par; nint=5
    ## 	 axp[0:1]=(1e+299,1e+308), usr[0:1]=(7.28752e+298,inf); i=9, ni=1
    pc <- 1:k
    matplot(np[iOk]/dFac, relE[iOk,], type=type, log=log, pch=pc, col=1+pc,
            main = "relative Errors  WRT  bd0(<mpfr-accurate>)",
            xlim = range(np)/dFac, # show full range
            xlab = paste0("np[iOk]", if(dFac != 1) sprintf("/ dFac,  dFac=%g",dFac)),
            ## could use  sfsmisc::pretty10exp(1e10, drop.1=TRUE)
            xaxt="n"); sfsmisc::eaxis(1, sub10=3)
    mtext(sprintf("bd0(x, np),  x = %g", x))
    if(k >= 2) legend("top", colnames(relE), pch=pc, lty=1:2, col=1+pc, bty="n")
    rug(np[!iOk]/dFac, col=2)
    axis(1, at=x/dFac, quote(x), col=2, col.axis=2, lwd=2, line=-1)
}

p.relE(bd0v.8)

## ==> FIXME:  a whole small (extreme) range where  bd0() is *better* than ebd0() !!!
with(bd0v.8, cbind(log2.lam = log2(np), np, relE)) ## around 2^[1018, 1021]


with(bd0v.8, stopifnot(exprs = {
    yhl[["yl"]] == 0 # which is not really good and should maybe change !
    ## Fixed now : both have 4 x Inf and then are equal {but do Note relE difference above!}
    all.equal(ebd0[["yh"]], bd0, tol = 4 * .Machine$double.eps)
}))
showProc.time()


## bd0()  and  ebd0()  are _Inf_  for smaller lambda's .. but they *must* be as true > DBL_MAX
## (almost: at the 4th value, Llam = 2^993, ideally the would *not* overflow: yM = 1.7598e+308)
bd0M <- bd0ver(x., Llam, mpfrPrec = 256, keepMpfr=TRUE)
with(bd0M, data.frame(log2.L = log2(np), bd0 = bd0, ebd0, bd0M. = format(bd0M, digits=8)))

matplot(log2(Llam), with(bd0M, cbind(bd0 = bd0/x., yh=ebd0[["yh"]]/x., asNumeric(bd0M/x.))),
        type="o", ylab = "*bd0*(x., L) / x.", pch=1:3,
        main= paste("ebd0(x., Lam) and bd0(*) for x=",format(x.)," *and* larg Lam"))
abline(h=0, lty=3, col=adjustcolor("gray20", 1/2))
axis(1, at=log2(x.), labels="log2(x.)", line=-1, tck=-1/16, col=2, col.axis=2)
legend("top", c("bd0()", "ebd0()", "MPFR bd0()"), bty="n", lty=1:3, pch=1:3, col=1:3)
dMax <- .Machine$double.xmax
abline(h = dMax / x., col=4, lty=3) ; ux <- par("usr")[1:2]
text(c(ux %*% c(3,1))/4, dMax/x., pos=3,
     sprintf("bd0(.) > DBL_MAX = %.5g", dMax), col=4)
showProc.time()


L.2 <- 2^c(seq(993, 1016, by=1/4), seq(1016+1/8, 1022, by=1/16))
bd0.2 <- bd0ver(x., L.2, mpfrPrec = 512, keepMpfr=TRUE)
p.relE(bd0.2)#, dFac=1)

if(!interactive()) # gets too expensive
    quit("no")

## zoom in more:
L.3 <- 2^c(seq(1019, 1021, by=1/128))
system.time(bd0.3 <- bd0ver(x., L.3, mpfrPrec = 1024, chkVerb=FALSE)) # 7.3 sec
p.relE(bd0.3) # up to 1e-11  rel.error !!

## different x :
system.time(bd0.2.2e307 <- bd0ver(2e307, L.2[L.2 > 1e306], mpfrPrec = 1024, chkVerb=FALSE)) # 1.2 s
p.relE(bd0.2.2e307)

## less large x .. still same problem:  ebd0() is worse than bd0()
L.4 <- 2^c(seq(1009, 1015, by=1/64))
system.time(bd0.2.2e305 <- bd0ver(2e305, L.4, mpfrPrec = 256, chkVerb=FALSE)) # 1.9 s
p.relE(bd0.2.2e305)

## less large x .. still same problem:  ebd0() is worse than bd0()
x <- 1e250
l2x <- round(64*log2(x))/64
str(L.x <- 2^(l2x + seq(-2,2, by=1/64)))
bd0.2.1e250 <- bd0ver(x, L.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0.2.1e250)

x <- 1e120
l2x <- round(256*log2(x))/256
str(L.x <- 2^(l2x + seq(-1,1, by=1/128)))
bd0.2.1e120 <- bd0ver(x, L.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0.2.1e120) # still rel.E -8e-13

x <- 1e20
l2x <- round(256*log2(x))/256
str(L.x <- 2^(l2x + seq(-1/2, 1/2, by=1/128)))
bd0.2.1e20 <- bd0ver(x, L.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0.2.1e20) #  rel.E -8e-12
apply(bd0.2.1e20$relE, 2, quantile)
##               ebd0         ebd0C           bd0
## 0%   -8.576108e-12 -8.576108e-12 -5.209563e-15
## 25%  -1.806152e-15 -1.806152e-15 -1.889226e-16
## 50%   3.757592e-16  3.757592e-16 -7.379701e-18
## 75%   4.833241e-15  4.833241e-15  1.476687e-16
## 100%  8.792719e-13  8.792719e-13  5.979845e-15

x <- 1e14
l2x <- round(256*log2(x))/256
str(L.x <- 2^(l2x + seq(-1/2, 1/2, by=1/128)))
bd0.2.1e14 <- bd0ver(x, L.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0.2.1e14) #  rel.E -2.7e-12
apply(bd0.2.1e14$relE, 2, quantile)
##               ebd0         ebd0C           bd0
## 0%   -2.675822e-12 -2.675822e-12 -4.087122e-15
## 25%  -2.227385e-15 -2.227385e-15 -2.372982e-16
## 50%   1.479060e-16  1.479060e-16  1.823222e-17
## 75%   2.234723e-15  2.234723e-15  2.340140e-16
## 100%  1.051075e-12  1.051075e-12  3.783036e-15


#      =
x <- 1e9
l2x <- round(256*log2(x))/256
str(L.x <- 2^(l2x + seq(-1/2, 1/2, by=1/128)))
bd0.2.1e9 <- bd0ver(x, L.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0.2.1e9) #  rel.E -2.7e-12
apply(bd0.2.1e9$relE, 2, quantile)
##               ebd0         ebd0C           bd0
## 0%   -5.279663e-12 -5.279663e-12 -5.611493e-15
## 25%  -1.377207e-15 -1.377207e-15 -1.696762e-16
## 50%  -2.192488e-16 -2.192488e-16  3.745953e-18
## 75%   1.578410e-15  1.578410e-15  1.491036e-16
## 100%  1.615619e-12  1.615619e-12  5.044319e-15

#      =
x <- 1e6
l2x <- round(256*log2(x))/256
str(L.x <- 2^(l2x + seq(-1, 1, by=1/128)))
bd0.2.1e6 <- bd0ver(x, L.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0.2.1e6) #  rel.E +- 4e-12
apply(bd0.2.1e6$relE, 2, quantile)
##               ebd0         ebd0C           bd0
## 0%   -4.034364e-13 -4.034364e-13 -4.774078e-15
## 25%  -2.305490e-15 -2.305490e-15 -1.934393e-16
## 50%  -3.393828e-17 -3.393828e-17 -1.192525e-17
## 75%   1.556812e-15  1.556812e-15  1.849302e-16
## 100%  4.098968e-12  4.098968e-12  4.140645e-15

## - {number of correct digits}:
matplot(L.x, log10(abs(bd0.2.1e6$relE[,-1])), type="l")


apply(bd0.2.1e6$relE, 2, quantile)


x <- 1e3
l2x <- round(256*log2(x))/256
str(L.x <- 2^(l2x + seq(-1/2, 1/2, by=1/128)))
bd0.2.1e3 <- bd0ver(x, L.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0.2.1e3) #  rel.E +- 4e-12
apply(bd0.2.1e3$relE, 2, quantile)
## - {number of correct digits}:
matplot(L.x, log10(abs(bd0.2.1e3$relE[,-1])), type="l")


##  now extend the "lambda" / np range:
str(np.x <- seq(.5, 2*x, length.out=1001))
bd0..1e3 <- bd0ver(x, np.x, mpfrPrec = 256, chkVerb=FALSE)
p.relE(bd0..1e3, log="") #  rel.E +- 4e-12
apply(bd0..1e3$relE, 2, quantile)
## even here, bd0(x, np)  is more accurate around  np ~= x
## but in the flanks,  ebd0 is better :
matplot(np.x, log10(abs(bd0..1e3$relE[,-1])), type="l")
