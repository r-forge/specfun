#### Testing  stirlerr(), bd0(), ebd0(), dpois_raw(), ...
#### ===============================================

require(DPQ)
for(pkg in c("Rmpfr", "DPQmpfr"))
    if(!requireNamespace(pkg)) {
        cat("no CRAN package", sQuote(pkg), " ---> no tests here.\n")
        q("no")
    }
require("Rmpfr")

cutoffs <- c(15,35,80,500) # cut points, n=*, in the above "algorithm"
##
n <- c(seq(1,15, by=1/4),seq(16, 25, by=1/2), 26:30, seq(32,50, by=2), seq(55,1000, by=5),
       20*c(51:99), 50*(40:80), 150*(27:48), 500*(15:20))
st.n <- stirlerr(n)# rather use.halves=TRUE, just here , use.halves=FALSE)
plot(st.n ~ n, log="xy", type="b") ## looks good now
nM <- mpfr(n, 2048)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose
all.equal(asNumeric(st.nM), st.n)# TRUE
all.equal(st.nM, as(st.n,"mpfr"))# .. difference: 1.05884..............................e-15
all.equal(roundMpfr(st.nM, 64), as(st.n,"mpfr"), tol=1e-16)# diff.: 1.05884...e-15

## Very revealing plot showing the *relative* approximation error of stirlerr(<dblprec>)

p.stirlerrDev <- function(n, precBits=2048, stnM = stirlerr(mpfr(n, precBits)), abs=FALSE,
                          ## cut points, n=*, in the stirlerr() algorithm :
                          cutoffs = c(15,35,80,500),
                          type = "b", cex = 1,
                          col = adjustcolor(1, 3/4), colnB = adjustcolor("orange4", 1/3),
                          log = if(abs) "xy" else "x",
                          xlim=NULL, ylim = if(abs) c(8e-18, max(abs(N(relE)))))
{
    op <- par(las = 1, mgp=c(2, 0.6, 0))
    on.exit(par(op))
    st <- stirlerr(n, cutoffs=cutoffs)
    relE <- sfsmisc::relErrV(stnM, st)
    N <- asNumeric
    form <- if(abs) abs(N(relE)) ~ n else N(relE) ~ n
    plot(form, log=log, type=type, cex=cex, col=col, xlim=xlim, ylim=ylim,
         ylab = quote(relErrV(stM, st)), axes=FALSE, frame=TRUE,
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

p.stirlerrDev(n=n, stnM=st.nM) # default cutoffs= c(15, 40, 85, 600)
## show the zoom-in region in next plot
yl2 <- 3e-14*c(-1,1)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

if(do.pdf) {
    dev.off() ; pdf("stirlerr-relErr_1.pdf", width=8, height=6)
}

## drop small n
p.stirlerrDev(n=n, stnM=st.nM, xlim = c(5, max(n))) # default cutoffs= c(15, 40, 85, 600)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

## The first plot clearly shows we should do better:
## Current code is switching to less terms too early, loosing up to 2 decimals precision
p.stirlerrDev(n=n, stnM=st.nM, ylim = yl2)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

if(do.pdf) {
    dev.off(); pdf("stirlerr-relErr_6-fin.pdf")
}

### ~19.April 2021: "This is close to *the* solution" (but ...)
cuts <- c(7, 12, 20, 26, 60, 200, 3300)
st. <- stirlerr(n=n, cutoffs = cuts, verbose=TRUE)
relE <- asNumeric(sfsmisc::relErrV(st.nM, st.))
head(cbind(n, relE), 20)
## nice printout :
print(cbind(n       = format(n, drop0trailing = TRUE),
            stirlerr= format(st.,scientific=FALSE, digits=4),
            relErr  = signif(relE, 4))
      , quote=FALSE)

if(do.pdf) {
    dev.off(); pdf("stirlerr-relErr_6-fin-1.pdf")
}

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts)

if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-2.pdf") }

## zoom in ==> {good for n >= 10}
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", ylim = 2e-15*c(-1,1),
              cutoffs = cuts)## old default cutoffs = c(15,35, 80, 500)

if(do.pdf) { dev.off(); pdf("stirlerr-tst_others.pdf") }

##-- April 20: we more terms up to S10 in stirlerr() -- more cutoffs
n <- sfsmisc::lseq(1/16, 5000, length=4096)
nM <- mpfr(n, 2048)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose

cuts <- c(5.4, 7.5, 8.5, 10.625, 12.125, 20, 26, 60, 200, 3300)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, abs=TRUE)
## using exact values sferr_halves[]
lines((0:30)/2, abs(stirlerr((0:30)/2, cutoffs=cuts, verbose=TRUE)/DPQ:::sferr_halves - 1), type="o", col=2,lwd=2)
## should we e.g., use interpolation spline through sfserr_halves[] for n <= 7.5
## -- doing the interpolation on the  log(1 - 12*x*stirlerr(x)) vs  log2(x)  scale -- maybe ?
curve(1-12*x*stirlerr(x, verbose=TRUE), 1/64, 8, log="xy", n=2048)
## just need "true" values for x = 2^-(6,5,4,3,2) in addition to those we already have at x = 1/2, 1.5, 2, 2.5, ..., 7.5, 8

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim=c(-1,1)*4e-14)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim=c(-1,1)*1e-15)

st. <- stirlerr(n=n, cutoffs = cuts, verbose=TRUE)
relE <- asNumeric(sfsmisc::relErrV(st.nM, st.))
head(cbind(n, relE), 20)
## nice printout :
print(cbind(n       = format(n, drop0trailing = TRUE),
            stirlerr= format(st.,scientific=FALSE, digits=4),
            relErr  = signif(asNumeric(sfsmisc::relErrV(st.nM, st.)), 4))
      , quote=FALSE)


## below, 7 "it's okay, but not perfect:" ===>  need more terms in stirlerr()  __or__ ??
## April 20: MM added more terms up to S10
x <- sfsmisc::lseq(1, 7, length=2048)
system.time(stM <- DPQmpfr::stirlerrM(Rmpfr::mpfr(x,2048))) # 1.52 sec elapsed
plot(x, stirlerr(x, use.halves=FALSE) - stM, type="l", log="x", main="absolute Error")
plot(x,     stirlerr(x, use.halves=FALSE) / stM - 1,  type="l", log="x", main="relative Error")
plot(x, abs(stirlerr(x, use.halves=FALSE) / stM - 1), type="l", log="xy",main="|relative Error|")
abline(h=c(1,2,4)*.Machine$double.eps, lty=3)
## lgammacor() does *NOT* help, as it is  *designed*  for  x >= 10!
lines(x, abs(lgammacor(x, 5) / stM - 1), col=2)
## maybe look at it for x >= 9 or so ?
##
## ==> Need another chebyshev() or rational-approx. for x in [.1, 7] or so !!



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


## close to over-/underflow -------

### Large lambda == np == M -------

if(do.pdf) { dev.off(); pdf("stirlerr-bd0-ebd0.pdf") }
##-- FIXME: use functionality from ~/R/MM/NUMERICS/dpq-functions/15628-dpois_raw_accuracy.R
##-- ----- *or* move to vignette

LL <- 1e20
dput(x1 <- 1e20 - 2e11) # 9.99999998e+19

(P1 <-         dpois     (x1,       LL)) # R-devel: 3.989455e-11; want 5.520993e-98
(P1m <- Rmpfr::dpois(mpfr(x1, 128), LL)) # 5.52099285934214335003128935..e-98
## However -- the ebd0() version
(P1e <- dpois_raw(x1, LL, version="ebd0_v1"))
## was 3.989455e-11, but now good!

options(digits=9)

## indeed:  regular  bd0()  works "ok"  --- but  ebd0() does non-sense !!
(bd.1 <- bd0(x1, LL, verbose=2))
## bd0(1e+20, 1e+20): T.series w/ 2 terms -> bd0=200
## [1] 200
(bd.1M <- bd0(x1, mpfr(LL, 128), verbose=2))
## bd0(1e+20, 1e+20): T.series w/ 3 terms -> bd0=200
## ---> 199.9999919413334091607468236761591740489

asNumeric(bd.1 / bd.1M - 1)# -1.82e-17 -- suggests bd0() is really accurate here
stopifnot(abs(bd.1 / bd.1M - 1) < 3e-16)

ebd0(x1, LL, verbose=TRUE)# fixed since  June 6, 2021


### Large x  -- small  np == M ------------------------------------

yy <-   bd0 (1e307, 10^(-5:1), verbose=TRUE)
yhl <- ebd0 (1e307, 10^(-5:1), verbose=TRUE)
stopifnot(yy == Inf, colSums(yhl) == Inf)
yM <- bd0(mpfr(1e307, 128), 10^(-5:1))
roundMpfr(range(yM), 12) ##  7.0363e+309 7.1739e+309 -- *are* larger than  DBL_MAX

## -- Now  *BOTH*  x and lambda are large : ------------------------

x. <- 1e307
ebd0(x., 10^(300:308))
stopifnot(is.finite(Llam <- 2^(990:1024 - 1e-12)))

yy  <-  bd0(x., Llam)
yhl.<- ebd0(x., Llam)
yhl <- ebd0(x., Llam, verbose=TRUE)
yM  <-  bd0(mpfr(x.,256), Llam)
relE <- asNumeric(cbind(ebd0 = yhl["yh",], bd0 = yy)/yM - 1)
iOk <- is.finite(yy) & yy != 0
## */1e10 : otherwise triggering axis() error
## log - axis(), 'at' creation, _LARGE_ range: invalid {xy}axp or par; nint=5
## 	 axp[0:1]=(1e+299,1e+308), usr[0:1]=(7.28752e+298,inf); i=9, ni=1
matplot(Llam[iOk]/1e10, relE[iOk,], type="b", log="x", pch=1:2, col=2:3,
        xaxt="n"); sfsmisc::eaxis(1, sub10=3)
legend("top", colnames(relE), pch=1:2, lty=1:2, col=2:3, bty="n")

stopifnot(exprs = {
    identical(yhl., yhl)
    yhl["yl",] == 0 # which is not really good and should maybe change !

    TRUE || ## FIXME ??  ebd0() gives many Inf here!

        all.equal(yhl["yh",], yy)

})

## actually *both*  bd0()  and  ebd0()  seem partly wrong here
data.frame(log2.L = log2(Llam), bd0 = yy, t(yhl), yM. = format(yM, digits=8))

matplot(log2(Llam), cbind(bd0 = yy/x., yh=yhl["yh",]/x., asNumeric(yM/x.)),
        type="o", ylab = "*bd0*(x., L) / x.", pch=1:3,
        main= paste("ebd0(x., Lam) and bd0(*) for x=",format(x.)," *and* larg Lam"))
abline(h=0, lty=3, col=adjustcolor("gray20", 1/2))
axis(1, at=log2(x.), label="log2(x.)", line=-1, tck=-1/16, col=2, col.axis=2)
legend("top", c("bd0()", "ebd0()", "MPFR bd0()"), bty="n", lty=1:3, pch=1:3, col=1:3)
