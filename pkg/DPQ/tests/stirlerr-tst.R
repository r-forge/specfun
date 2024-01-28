#### Testing  stirlerr()
#### =================== {previous 2nd part of this, now -->>> ./bd0-tst.R <<<---
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

do.pdf <- TRUE # (manually)
do.pdf <- !dev.interactive(orNone = TRUE)
do.pdf
if(do.pdf) {
    pdf.options(width = 8, height = 6) # for all pdf plots
    pdf("stirlerr-tst.pdf")
}

showProc.time()
(doExtras <- DPQ:::doExtras())
(noLdbl <- (.Machine$sizeof.longdouble <= 8)) ## TRUE when --disable-long-double
options(width = 100, nwarnings = 1e5)

abs19 <- function(r) pmax(abs(r), 1e-19) # cut  |err| to positive {for log-plots}

## to enhance  |rel.Err| plots:  {also in ./pow-tst.R and  ~/R/Pkgs/Rmpfr/tests/special-fun-ex.R }
drawEps.h <- function(p2 = -(53:51), negative=FALSE, side = 4, lty = 3, lwd = 2, col = adjustcolor(2, 1/2)) {
    twop <- if(negative) c(outer(2^p2, c(-1,1))) else 2^p2
    labL <- lapply(p2, function(p) substitute(2^E, list(E=p)))
    if(negative) labL <- c(lapply(p2, function(p) substitute(-2^E, list(E=p))), labL)
    abline(h = twop, lty=lty, lwd=lwd, col=col)
    if(negative) abline(h=0, col=adjustcolor("gray20", 1/2), lwd=max(1, lwd))# at least thick as others
    axis(side, las=2, line=-1, at = twop, labels = as.expression(labL),
         col.axis = col, col=NA, col.ticks=NA)
}

cutoffs <- c(15,35,80,500) # cut points, n=*, in the stirlerr() "algorithm"
##
n <- c(seq(1,15, by=1/4),seq(16, 25, by=1/2), 26:30, seq(32,50, by=2), seq(55,1000, by=5),
       20*c(51:99), 50*(40:80), 150*(27:48), 500*(15:20))
st.n <- stirlerr(n, "R3")# rather use.halves=TRUE; but here , use.halves=FALSE
plot(st.n ~ n, log="xy", type="b") ## looks good now (straight line descending {in log-log !}
nM <- mpfr(n, 2048)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose
all.equal(asNumeric(st.nM), st.n)# TRUE
all.equal(st.nM, as(st.n,"mpfr"))# .. difference: 3.381400e-14 was 1.05884.........e-15
all.equal(roundMpfr(st.nM, 64), as(st.n,"mpfr"), tolerance=1e-16)# (ditto)


## --- Look at the direct formula -- why is it not good for n ~= 5 ?
##
## Preliminary Conclusions :
## 1. there is *some* cancellation even for small n (how much?)
## 2. lgamma1p(n) does really not help much compared to lgamma(n+1) --- but a tiny bit in some cases

###  1. Investigating  lgamma1p(n)  vs  lgamma(n+1)  for n < 1 =============================================

##' @title Relative Error of lgamma(n+1) vs lgamma1p() vs MM's stirlerrD2():
##' @param n numeric, typically n << 1
##' @param precBits
##' @return relative error WRT  mpfr(n, precBits)
##' @author Martin Maechler
relE.lgam1 <-  function(n, precBits = if(doExtras) 1024 else 320) {
    M_LN2PI <- 1.837877066409345483 # ~ log(2*Const("pi",60)); very slightly more accurate than log(2*pi)
    st  <- lgamma(n +1)  - (n +0.5)*log(n)  + n  - M_LN2PI/2
    st. <- lgamma1p(n)   - (n +0.5)*log(n)  + n  - M_LN2PI/2      # "lgamma1p"
    st2 <- lgamma(n) + n*(1-(l.n <- log(n))) + (l.n - M_LN2PI)/2  # "MM2"
    st0 <- -(l.n + M_LN2PI)/2                                     # "n0"
    nM <- mpfr(n, precBits)
    stM <- lgamma(nM+1)  - (nM+0.5)*log(nM) + nM - log(2*Const("pi", precBits))/2
    ## stM <- roundMpfr(stM, 128)
    cbind("R3"      = asNumeric(relErrV(stM, st))
        , "lgamma1p"= asNumeric(relErrV(stM, st.))
        , "MM2"	    = asNumeric(relErrV(stM, st2))
        , "n0"      = asNumeric(relErrV(stM, st0))
          )
}

n <- 2^-seq.int(1022, 1, by = -1/4)
relEx <- relE.lgam1(n)

## Is *equivalent* to 'new' stirlerr_simpl(n version = *):
stir.allS <- function(n) sapply(eval(formals(stirlerr_simpl)$version),
                                function(v) stirlerr_simpl(n, version=v))
stirS <- stir.allS(n)
releS <- asNumeric(relErrV(stirlerr(mpfr(n,  256)), stirS)) # prec = 256 should suffice!
all.equal(relEx, releS, tolerance = 0) # see TRUE on Linux
stopifnot(all.equal(relEx, releS, tolerance = 1e-15))

cols <- c("gray30", adjustcolor(c(2,3,4), 1/2));  lwd <- c(1, 3,3,3)
stopifnot((k <- length(cols)) == ncol(relEx), k == length(lwd))
matplot(n, relEx, type = "l", log="x", col=cols, lwd=lwd, ylim = c(-1,1)*4.5e-16,
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
mtext("really small errors are dominated by small (< 2^-53) errors of log(n)")
## very interesting: there are different intervals  <---> log(n) Qpattern !!
## -- but very small difference, only for n >~= 1/1000  but not before
drawEps.h(negative=TRUE) # abline(h= c(-4,-2:2, 4)*2^-53, lty=c(2,2,2, 1, 2,2,2), col="gray")
legend("topleft", legend = colnames(relEx), col=cols, lwd=3)

## zoomed in a bit:
n. <- 2^-seq.int(400,2, by = -1/4)
relEx. <- relE.lgam1(n.)
matplot(n., relEx., type = "l", log="x", col=cols, lwd=lwd, ylim = c(-1,1)*4.5e-16,
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
drawEps.h(negative=TRUE)
legend("topleft", legend = colnames(relEx.), col=cols, lwd=3)

##====> Absolute errors (and look at  "n0") --------------------------------------
matplot(n., abs19(relEx.), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(4e-17, 5e-16),
        main = quote(abs(relErr(stirlerr_simpl(n, '*')))))
drawEps.h(); legend("top", legend = colnames(relEx.), col=cols, lwd=3)
lines(n., abs19(relEx.[,"n0"]), type = "o", cex=1/4, col=cols[4], lwd=2)

## more zooom-in
n.2 <- 2^-seq.int(85, 50, by= -1/100)
stirS.2 <- sapply(c("R3", "lgamma1p", "n0"), function(v) stirlerr_simpl(n.2, version=v))
releS.2 <- asNumeric(relErrV(stirlerr(mpfr(n.2, 320)), stirS.2))

matplot(n.2, abs19(releS.2), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(4e-17, 5e-16),
        main = quote(abs(relErr(stirlerr_simpl(n, '*')))))
drawEps.h(); legend("top", legend = colnames(releS.2), col=cols, lwd=3)
abline(v = 5e-17, col=(cb <- adjustcolor("skyblue4", 1/2)), lwd=2, lty=3)
axis(1, at=5e-17, col.axis=cb, line=-1/4, cex = 3/4)

matplot(n.2, abs19(releS.2), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(4e-17, 5e-16),
        xaxt="n", xlim = c(8e-18, 1e-15), ## <<<<<<<<<<<<<<<<<<< Zoom-in
        xlab = quote(n), main = quote(abs(relErr(stirlerr_simpl(n, '*')))))
eaxis(1); drawEps.h(); legend("top", legend = colnames(releS.2), col=cols, lwd=3)
abline(v = 5e-17, col=(cb <- adjustcolor("skyblue4", 1/2)), lwd=2, lty=3)
mtext('stirlerr_simpl(*, "n0")  is as good as others for n <= 5e-17', col=adjustcolor(cols[3], 2))
##     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ===> "all *but* "n0" approximations for "larger" small n:
n2 <- 2^-seq.int(20,0, length.out=1000)
relEx2 <- relE.lgam1(n2)[,c("R3", "lgamma1p", "MM2")] # "n0" is "bad": |relE| >= 2.2e-6 !
cols <- c("gray30", adjustcolor(c(2,4), 1/2));  lwd <- c(1, 3,3)
stopifnot((k <- length(cols)) == ncol(relEx2), k == length(lwd))
matplot(n2, relEx2, type = "l", log="x", col=cols, lwd=lwd, ylim = c(-3,3)*1e-15, xaxt="n",
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
eaxis(1, sub10=c(-3,0)); drawEps.h(negative=TRUE)
legend("topleft", legend = colnames(relEx2), col=cols, lwd=3)
##==> "MM" is *NOT* good for n < 1  *but*

## "the same" -- even larger small n:
n3 <- seq(.01, 5, length=1000)
relEx3 <- relE.lgam1(n3)[,c("R3", "lgamma1p", "MM2")] # "no" is "bad" ..
stopifnot((k <- length(cols)) == ncol(relEx3), k == length(lwd))

matplot(n3, relEx3, type = "l", col=cols, lwd=lwd,
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
legend("topleft", legend = colnames(relEx3), col=cols, lwd=3)
drawEps.h(negative=TRUE)

matplot(n3, abs19(relEx3), type = "l", col=cols, lwd=lwd,
        log="y", ylim = 2^-c(54, 44), yaxt = "n", ylab = quote(abs(relE)), xlab=quote(n),
        main = "|relative errors| of direct (approx.) formula for stirlerr(n), small n")
eaxis(2, cex.axis=0.9); legend("topleft", legend = colnames(relEx3), col=cols, lwd=3)
drawEps.h()
lines(n3, smooth.spline(abs(relEx3)[,1], df=12)$y, lwd=3, col=cols[1])
lines(n3, smooth.spline(abs(relEx3)[,2], df=12)$y, lwd=3, col=adjustcolor(cols[2], 1/2))
lines(n3, smooth.spline(abs(relEx3)[,3], df=12)$y, lwd=4, col=adjustcolor(cols[3], offset = rep(.2,4)))
## ===>  from  n >~= 1,  "MM2" is definitely better up to n = 5 !!

## Check  log() only :
plot(n, asNumeric(relErrV(log(mpfr(n, 256)), log(n))), ylim = c(-1,1)*2^-53,
     log="x", type="l", xaxt="n") ## ===> indeed --- log(n) approximation pattern !!
eaxis(1) ; drawEps.h(negative=TRUE)

## =========== "R3"  vs  "lgamma1p" -------------------------- which is better?

## really for the very small n, all is dominated by -(n+0.5)*log(n);  and lgamma1p() is unnecessary!
i <- 1:20; ni <- n[i]
lgamma1p(ni)
- (ni +0.5)*log(ni)  + ni

## much less extreme:
n2 <- lseq(2^-12, 1/2, length=1000)
relE2 <- relE.lgam1(n2)[,-4]

cols <- c("gray30", adjustcolor(2:3, 1/2)); lwd <- c(1,3,3)
matplot(n2, relE2, type = "l", log="x", col=cols, lwd=lwd)
legend("topleft", legend=colnames(relE2), col=cols, lwd=2, lty=1:3)
drawEps.h(negative=TRUE)

matplot(n2, abs19(relE2), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(6e-17, 1e-15),
        xaxt = "n"); eaxis(1, sub10=c(-2,0))
legend("topleft", legend=colnames(relE2), col=cols, lwd=2, lty=1:3)
drawEps.h()
## "MM2" is *worse* here, n < 1/2
for(j in 1:3) lines(n2, smooth.spline(abs(relE2[,j]), df=10)$y, lwd=3,
                    col=adjustcolor(cols[j], 1.5, offset = rep(-1/4, 4)))
## "lgammap  very slightly better in [0.002, 0.05] ...
## "TODO": draw 0.90-quantile curves {--> cobs::cobs() ?} instead of mean-curves?

## which is better?  ... "random difference"
d.absrelE <- abs(relE2[,"R3"]) - abs(relE2[,"lgamma1p"])
plot (n2, d.absrelE, type = "l", log="x", # no clear picture ...
      main = "|relE_R3|  -  |relE_lgamma1p|", axes=FALSE, frame.plot=TRUE)
eaxis(1, sub10=c(-2,1)); eaxis(2); axis(3, at=max(n2)); abline(v = max(n2), lty=3, col="gray")
## 'lgamma1p' very slightly better:
lines(n2, smooth.spline(d.absrelE, df=12)$y, lwd=3, col=2)

## not really small n at all == here see, how "bad" the direct formula gets for 1 < n < 10 or so
n3 <- lseq(2^-14, 2^2, length=800)
relE3 <- relE.lgam1(n3)[, -4]

matplot(n3, relE3, type = "l", log="x", col=cols, lty=1, lwd = c(1,3),
        main = quote(rel.lgam1(n)), xlab=quote(n))

matplot(n3, abs19(relE3),  type = "l", log="xy", col=cols, lwd = c(1,3), xaxt="n",
        main = quote(abs(rel.lgam1(n))), xlab=quote(n), ylim = c(2e-17, 4e-14))
drawEps.h(); eaxis(1, sub10=c(-2,3))
legend("topleft", legend=colnames(relE3), col=cols, lwd=2)
## very small difference --- draw the 3 smoothers :
for(j in 1:3) {
    ll <- lowess(log(n3), abs19(relE3[, j]), f= 1/12)
    with(ll, lines(exp(x), y, col=adjustcolor(cols[j], 1.5), lwd=3))
}
## ==> lgamma1p(.) very slightly in n ~ 10^-4 -- 10^-2 --- but not where it matters: n ~ 0.1 -- 1 !!
##     "MM2"  gets best from  n >~ 1 !
abline(v=1, lty=3, col = adjustcolor(1, 3/4))


###  2. relErr( stirlerr(.) ) ============================================================

##' Very revealing plot showing the *relative* approximation error of stirlerr(<dblprec>)
##'
p.stirlerrDev <- function(n, precBits = if(doExtras) 2048L else 512L,
                          stnM = stirlerr(mpfr(n, precBits), use.halves=use.halves, verbose=verbose),
                          abs = FALSE,
                          ## cut points, n=*, in the stirlerr() algorithm; "FIXME": sync with ../R/dgamma.R <<<<
                          scheme = c("R3", "R4.4_0"),
                          cutoffs = switch(match.arg(scheme)
                                         , R3   = c(15, 35, 80, 500)
                                         , R4.4_0 = c(5.0, 5.3, 5.4, 5.7, 6.1, 6.5, 7.0, 7.9, 8.75, 10.5, 13,
                                                      20, 26, 60, 200, 3300, 17.4e6)
                                           ## {FIXME: need to sync} <==> ../man/stirlerr.Rd <==> ../R/dgamma.R
                                           ),
                          use.halves = missing(cutoffs),
                          direct.ver = c("R3", "lgamma1p", "MM2", "n0"),
                          verbose = getOption("verbose"),
                          type = "b", cex = 1,
                          col = adjustcolor(1, 3/4), colnB = adjustcolor("orange4", 1/3),
                          log = if(abs) "xy" else "x",
                          xlim=NULL, ylim = if(abs) c(8e-18, max(abs(N(relE)))))
{
    op <- par(las = 1, mgp=c(2, 0.6, 0))
    on.exit(par(op))
    require("Rmpfr"); require("sfsmisc")
    st <- stirlerr(n, scheme=scheme, cutoffs=cutoffs, use.halves=use.halves, direct.ver=direct.ver,
                   verbose=verbose)
    relE <- relErrV(stnM, st) # eps0 = .Machine$double.xmin
    N <- asNumeric
    form <- if(abs) abs(N(relE)) ~ n else N(relE) ~ n
    plot(form, log=log, type=type, cex=cex, col=col, xlim=xlim, ylim=ylim,
         ylab = quote(relErrV(stM, st)), axes=FALSE, frame.plot=TRUE,
         main = sprintf("stirlerr(n, cutoffs) rel.error [wrt stirlerr(Rmpfr::mpfr(n, %d))]",
                        precBits))
    eaxis(1, sub10=3)
    eaxis(2)
    mtext(paste("cutoffs =", deparse1(cutoffs)))
    ylog <- par("ylog")
    ## FIXME:  improve this ---> drawEps.h() above
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

showProc.time()

n <- lseq(2^-10, 1e10, length=4096)
n <- lseq(2^-10, 5000, length=4096)
## store "expensive" stirlerr() result, and re-use many times below:
nM <- mpfr(n, if(doExtras) 2048 else 512)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose

p.stirlerrDev(n=n, stnM=st.nM, use.halves = FALSE) # default cutoffs= c(15, 40, 85, 600)
p.stirlerrDev(n=n, stnM=st.nM, use.halves = FALSE, ylim = c(-1,1)*1e-12) # default cutoffs= c(15, 40, 85, 600)

## show the zoom-in region in next plot
yl2 <- 3e-14*c(-1,1)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

if(do.pdf) { dev.off() ; pdf("stirlerr-relErr_1.pdf") }

## drop n < 5:
p.stirlerrDev(n=n, stnM=st.nM, xlim = c(7, max(n)), use.halves=FALSE) # default cutoffs= c(15, 40, 85, 600)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

## The first plot clearly shows we should do better:
## Current code is switching to less terms too early, loosing up to 2 decimals precision
if(FALSE) # no visible difference {use.halves = T / F }:
p.stirlerrDev(n=n, stnM=st.nM, ylim = yl2, use.halves = FALSE)
p.stirlerrDev(n=n, stnM=st.nM, ylim = yl2, use.halves = TRUE)# exact at n/2 (n <= ..)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

showProc.time()


if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-1.pdf") }

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
## and zoom in:
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim = yl2)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim = yl2/20)

if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-2.pdf") }

## zoom in ==> {good for n >= 10}
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", ylim = 2e-15*c(-1,1),
              cutoffs = cuts)## old default cutoffs = c(15,35, 80, 500)

if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-3.pdf") }


##-- April 20: have more terms up to S10 in stirlerr() --> can use more cutoffs
n <- n5m <- lseq(1/64, 5000, length=4096)
nM <- mpfr(n, if(doExtras) 2048L # a *lot* accuracy for stirlerr(nM,*)
	      else 512L)
ct10.1 <- c(            5.4, 7.5, 8.5, 10.625, 12.125, 20, 26, 60, 200, 3300)# till 2024-01-19
ct10.2 <- c(            5.4, 7.9, 8.75,10.5  , 13,     20, 26, 60, 200, 3300)
cuts <-
ct12.1 <- c(5.22, 6.5,  7.0, 7.9, 8.75,10.5  , 13,     20, 26, 60, 200, 3300)
##        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5.25 is "too small" but the direct formula is already really bad there, ...
st.nM <- roundMpfr(stirlerr(nM, use.halves=FALSE, ## << on purpose;
                            verbose=TRUE), precBits = 128)
## NB: for x=xM <mpfr>; `cutoffs` are *not* used.
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0")
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0", ylim = c(-1,1)*3e-15)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0", ylim = c(-1,1)*1e-15)

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0", abs=TRUE)
axis(1,at= 2:6, col=NA, col.axis=(cola <- "lightblue"), line=-3/4)
abline(v = 2:6, lty=3, col=cola)
if(FALSE)## using exact values sferr_halves[] *instead* of MPFR ones: ==> confirmation they lay on top
lines((0:30)/2, abs(stirlerr((0:30)/2, cutoffs=cuts, verbose=TRUE)/DPQ:::sferr_halves - 1), type="o", col=2,lwd=2)

if(FALSE) ## nice (but unneeded) printout :
print(cbind(n       = format(n, drop0trailing = TRUE),
            stirlerr= format(st.,scientific=FALSE, digits=4),
            relErr  = signif(relE, 4))
      , quote=FALSE)

showProc.time()

## ========== where should the cutoffs be ? ===================================================

.stirl.cutoffs <- function(scheme)
    eval(do.call(substitute, list(formals(stirlerr)$cutoffs, list(scheme = scheme))))
drawCuts <- function(scheme, axis=NA, lty = 3, col = "skyblue", ...) {
    abline(v = (ct <- .stirl.cutoffs(scheme)), lty=lty, col=col, ...)
    if(is.finite(axis)) axis(axis, at = ct, col=col)
}
axisCuts <- function(scheme, side = 4, col = "skyblue", ...)
    axis(side, at = .stirl.cutoffs(scheme), col.axis = col, col=NA, col.ticks=NA, ...)
mtextCuts <- function(cutoffs, scheme, ...) {
    if(!missing(scheme)) cutoffs <- .stirl.cutoffs(scheme)
    mtext(paste("cutoffs =", deparse1(cutoffs)), ...)
}


if(do.pdf) { dev.off(); pdf("stirlerr-tst_order_k.pdf") }

mK <- 17L # := max(k)
## order = k = 1:mK  terms in series approx:
k <- 1:mK
n <- 2^seq(1, 28, by=1/16)
nM <- mpfr(n, 1024)
stnM <- stirlerr(nM) # the "true" values

stirlOrd <- sapply(k, function(k) stirlerr(n, order = k))
relE <- asNumeric(stirlOrd/stnM -1) # "true" relativ error

## use a "smooth" but well visible polette :
palROBG <- colorRampPalette(c("red", "darkorange2", "blue", "seagreen"), space = "Lab")
palette(adjustcolor(palROBG(mK), 3/4))

(tit.k  <- substitute(list(    stirlerr(n, order=k) ~~"error",  k == 1:mK),  list(mK = mK)))
(tit.kA <- substitute(list(abs(stirlerr(n, order=k) ~~"error"), k == 1:mK),  list(mK = mK)))

matplotB(n, relE, cex=2/3, ylim = c(-1,1)*1e-13, col=k,
        log = "x", xaxt="n", main = tit.k)
eaxis(1, nintLog = 20)
drawCuts("R4.4_0")

## zoom in (ylim)
matplotB(n, relE, cex=2/3, ylim = c(-1,1)*5e-15, col=k,
        log = "x", xaxt="n", main = tit.k)
eaxis(1, nintLog = 20); abline(h = (-2:2)*2^-53, lty=3, lwd=1/2)
drawCuts("R4.4_0")

## log-log  |rel.Err|  -- "linear"
matplotB(n, abs19(relE), cex=2/3, col=k, ylim = c(8e-17, 1e-3), log = "xy", main=tit.kA)
mtext(paste("k =", deparse(k))) ; abline(h = 2^-(53:51), lty=3, lwd=1/2)
drawCuts("R4.4_0")

## zoom into the critical n region
nc <- seq(3, 12, by=1/32)
ncM <- mpfr(nc, 1024)
stncM <- stirlerr(ncM) # the "true" values
stirlO.c <- sapply(k, function(k) stirlerr(nc, order = k))
relEc <- asNumeric(stirlO.c/stncM -1) # "true" relativ error

## log-log  |rel.Err|  -- "linear"
matplotB(nc, abs19(relEc), cex=2/3, col=k, ylim = c(2e-17, 1e-8),
        log = "xy", xlab = quote(n), main = quote(abs(relErr(stirlerr(n, order==k)))))
mtext(paste("k =", deparse(k))) ; abline(h = 2^-(53:51), lty=3, lwd=1/2)
lines(nc, abs19(asNumeric(stirlerr_simpl(nc)/stncM - 1)), lwd=3, col=adjustcolor(2, 2/3))
legend(10^par("usr")[1], 1e-9, legend=paste0("k=", k), bty="n", lwd=2,
       col=k, lty=1:5, pch= c(1L:9L, 0L, letters)[seq_along(k)])
drawCuts("R4.4_0")

##--- Accuracy of "R4.4_0" -------------------------------------------------------

nc <- seq(5, 200, by=1/256) # for a bigger pix
nc <- seq(4.75, 7, by=1/1024)
ncM <- mpfr(nc, 1024)
stncM <- stirlerr(ncM) # the "true" values
stirl.440 <- stirlerr(nc, scheme = "R4.4_0")
stirl.3   <- stirlerr(nc, scheme = "R3")
relE440 <- asNumeric(relErrV(stncM, stirl.440))
relE3   <- asNumeric(relErrV(stncM, stirl.3  ))

plot(nc, abs19(relE440), xlab=quote(n), main = quote(abs(relErr(stirlerr(n, '"R4.4_0"')))),
     type = "l", log = "xy", ylim = c(4e-17, 1e-13))
mtextCuts(scheme="R4.4_0")
drawCuts("R4.4_0", lty=2, lwd=2, axis=4)
drawEps.h()
if(max(nc) <= 10) abline(v = 5+(0:20)/10, lty=3, col=adjustcolor(4, 1/2))
if(TRUE) { # but just so ...
    c3 <- adjustcolor("royalblue", 1/2)
    lines(nc, pmax(abs(relE3), 1e-18), col=c3)
    title(quote(abs(relErr(stirlerr(n, '"R3"')))), adj=1, col.main = c3)
    drawCuts("R3", lty=4, col=c3); mtextCuts(scheme="R3", adj=1, col=c3)
}
addOrd <- TRUE
addOrd <- FALSE
if(addOrd) {
    i <- (15 <= nc & nc <= 85)
    ni <- nc[i]
    for(k in 7:17) lines(ni, abs19(asNumeric(relErrV(stncM[i], stirlerr(ni, order=k)))), col=k)
}

## ------ stirlerr(. order = *)      [again? -- keep?]
nc <- seq(4.75, 7, by=1/1024)
ncM <- mpfr(nc, 1024)
stncM <- stirlerr(ncM) # the "true" values
stirlO.c <- sapply(k, function(k) stirlerr(nc, order = k))
relEc <- asNumeric(stirlO.c/stncM -1) # "true" relativ error

## log  |rel.Err|  -- "linear"
matplotB(nc, pmax(abs(relEc), 1e-19), col=k, cex = 2/3, ylim = c(2e-17, 1e-11),
        log = "y", xlab = quote(n), main = quote(abs(relErr(stirlerr(n, order==k)))))
mtext(paste("k =", deparse(k))) ; abline(h = 2^-(53:51), lty=3, lwd=1/2)
lines(nc, pmax(asNumeric(stirlerr_simpl(nc)/stncM - 1), 1e-19), lwd=1.5, col=adjustcolor(2, 0.4))
k. <- k[k >= 6]
legend("bottomleft", legend=paste0("k=", k.), bty="n", lwd=2,
       col=k., lty=1:5, pch= c(1L:9L, 0L, letters)[k.])
abline(v = 5+(0:20)/10, lty=3, col=adjustcolor(10, 1/2))
drawCuts("R4.4_0")


## ========== Try a slightly better direct formula ======================================================

## after some trial error:  gamma(n+1) = n*gamma(n) ==> lgamma(n+1) = lgamma(n) + log(n)
## stirlerrD2 <- function(n) lgamma(n) + n*(1-(l.n <- log(n))) + (l.n - log(2*pi))/2

## the above plot
lines(nc, pmax(asNumeric(stirlerr_simpl(nc, "MM2")/stncM - 1), 1e-19), lwd=1, col=adjustcolor(10, 0.4))

palette("default")

if(do.pdf) { dev.off(); pdf("stirlerr-tst_others.pdf") }

n <- n5m # (1/64 ... 5000) -- goes with st.nM
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, abs=TRUE)
axis(1,at= 2:6, col=NA, col.axis=(cola <- "lightblue"), line=-3/4)
abline(v = 2:6, lty=3, col=cola)
i.n <- 1 <= n & n <= 15 ; cr2 <- adjustcolor(2, 0.6)
lines(n[i.n], abs(relErrV(st.nM[i.n], stirlerr_simpl(n[i.n], "MM2" ))), col=cr2, lwd=2)
legend(20, 1e-13, legend = quote(  relErr(stirlerr_simpl(n, '"MM2"'))), col=cr2, lwd=3, bty="n")


n <- seq(1,6, by= 1/200)
stM <- stirlerr(mpfr(n, 512), use.halves = FALSE)
relE <- asNumeric(relErrV(stM, cbind(stD = stirlerr_simpl(n), stD2 = stirlerr_simpl(n, "MM2"))))
signif(apply(abs(relE), 2, summary), 4)
##               stD      stD2
## Min.    3.093e-17 7.081e-18
## 1st Qu. 2.777e-15 1.936e-15
## Median  8.735e-15 5.238e-15
## Mean    1.670e-14 1.017e-14  .. well "67% better"
## 3rd Qu. 2.380e-14 1.408e-14
## Max.    1.284e-13 9.382e-14

c2 <- adjustcolor(1:2, 0.6)
matplotB(n, pmax(abs(relE), 1e-19), type="o", cex=3/4, log="y", ylim = c(8e-17, 1.3e-13), yaxt="n", col=c2)
eaxis(2)
abline(h = 2^-53, lty=1, col="gray")
smrelE <- apply(abs(relE), 2, \(y) lowess(n, y, f = 0.1)$y)
matlines(n, smrelE, lwd=3, lty=1)
legend("topleft", legend = expression(stirlerr_simpl(n), stirlerr_simple(n, "MM2")),
       bty='n', col=c2, lwd=3, lty=1)
drawEps.h(-(53:48))


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




##=== Really, dpois_raw()  and  dbinom_raw()  *both* use stirlerr(x)  for "all"  'x > 0'
##            ~~~~~~~~~~~       ~~~~~~~~~~~~             ===========              ===== !

## below, 6 "it's okay, but *far* from perfect:" ===>  need more terms in stirlerr() [
## April 20: MM added more terms up to S10; 2024-01: up to S12 ..helps a little only
x <- lseq(1/16, 6, length=2048)
system.time(stM <- DPQmpfr::stirlerrM(Rmpfr::mpfr(x,2048))) # 1.7 sec elapsed
plot(x, stirlerr(x, use.halves=FALSE) - stM,      type="l", log="x", main="absolute Error")
plot(x, stirlerr(x, use.halves=FALSE) / stM - 1,  type="l", log="x", main="relative Error")
plot(x, abs(stirlerr(x, use.halves=FALSE) / stM - 1), type="l", log="xy",main="|relative Error|")
abline(h=c(1,2,4)*.Machine$double.eps, lty=3)
## lgammacor() does *NOT* help, as it is  *designed*  for  x >= 10!
##
## ==> Need another chebyshev() or rational-approx. for x in [.1, 7] or so !!

##=============> For now, see ../Misc/stirlerr-trms.R  <===============
##                            ~~~~~~~~~~~~~~~~~~~~~~~
showProc.time()
