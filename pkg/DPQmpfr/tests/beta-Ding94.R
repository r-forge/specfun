require(DPQmpfr)
library ( Rmpfr)

options(warn = 1)# warnings *immediately*
(doExtras <- DPQmpfr:::doExtras())

##  *.time()  utilities   from
if(FALSE)
source(system.file("test-tools-1.R", package="Matrix"), keep.source=FALSE)
##    MM = ~/R/Pkgs/Matrix/inst/test-tools-1.R

showSys.time <- function(expr, ...) {
    ## prepend 'Time' for R CMD Rdiff
    st <- system.time(expr, ...)
    writeLines(paste("Time", capture.output(print(st))))
    invisible(st)
}
showProc.time <- local({ ## function + 'pct' variable
    pct <- proc.time()
    function(final="\n") { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	## 'Time ..' *not* to be translated:  tools::Rdiff() skips its lines!
	cat('Time elapsed: ', (pct - ot)[1:3], final)
    }
})


### pbeta*() :

pbetaD94(mpfr(0.864,  140),  5, 5, 54, eps=1e-10)# n=233 at convergence, 0.456302619295....
pbetaD94(mpfr(0.864,  140),  5, 5, 54, eps=1e-30)# n=572  "      "     , 0.45630261933697895485..

pbetaD94(mpfr(0.956,  140),  5, 5, 170, eps=1e-10)# n= 765, 0.602242164945
pbetaD94(mpfr(0.956,  140),  5, 5, 170, eps=1e-20)# n=1323, 0.60224216500116548

pbetaD94(mpfr(0.8686, 140), 10,10,  54, eps=1e-10)# n=304, 0.918779110843..
pbetaD94(mpfr(0.8686, 140), 10,10,  54, eps=1e-20)# n=497, 0.9187791109260769059699

pbetaD94(mpfr(0.8787, 140), 20,20,  54, eps=1e-10)# n=456, 0.99986765729631163
pbetaD94(mpfr(0.8787, 180), 20,20,  54, eps=1e-20)# n=

pbetaD94(mpfr(0.922,  140), 20,20, 250, eps=1e-10)# n=743, 0.964119072835964766
pbetaD94(mpfr(0.922,  180), 20,20, 250, eps=1e-20)# n=

## MM:  Try large lambdas ==> what we need is a *RELATIVE* bound  --> now implemented
pbetaD94(mpfr(0.93,   160), 20,25, 10000)            # n= 5175 1.008067955986..e-116
pbetaD94(mpfr(0.93,   160), 20,25, 10000, eps=1e-20) # n= 5513 1.008067956082285281372...e-116
pbetaD94(mpfr(0.93,   256), 20,25, 10000, eps=1e-30) # n= 5850 1.008067956082285281381936816362733357442687434229...e-116
## double the precBits only
pbetaD94(mpfr(0.93,   512), 20,25, 10000, eps=1e-30) # n= 5850 1.008067956082285281381936816362733357442687434229...e-116
## double the precBits again + eps !!
pbetaD94(mpfr(0.93,  1024), 20,25, 10000, eps=1e-50) # n= 6520 1.008067956082285281381936816363684457888595749378086614..e-116

## and much more extreme accuracy.. seems to work, too:
pbetaD94(mpfr(0.93,  1024), 20,25, 10000, eps=1e-150)# n= 9828 1.00806795608228528138193681636368445788859574937809624269896.......e-116


### dbeta*() :

" TODO "

### qbeta*() --- careful about speed / time usage

## This now works, thanks to log_scale !
(q.5.250 <- qbetaD94(0.95, .5, 250, ncp = 0))

## but using mpfr of course works:
require("Rmpfr") # (is in DPQmpfr's strict dependencies)
qbetaD94(1 - 1/mpfr(20 ,64), .5, 250)                           # 0.00766110246787361539366
qbetaD94(1 - 1/mpfr(20,128), .5, 250, eps = 1e-30, delta=1e-20) # 0.0076611024678614272722717848...

## Compare with  Table 3  of  Baharev_et_al 2017 %% ===> ../man/qbBaha2017.Rd <<<<<<<<<<<<
aa <- c(0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 25)
bb <- c(1:15, 10*c(2:5, 10, 25, 50))

## Choose one of these
p <- 0.95 # = 1 - alpha <<<--- using double precision everywhere below
p <- 1 - 1/mpfr(20,128) # = 1 - alpha <<<--- using MPFR everywhere below

delta <- 1e-18
eps <- delta^2 # default; (delta,eps) are needed for safe-file name
qbet <- matrix(NA_real_, length(aa), length(bb),
               dimnames = list(a = formatC(aa), b = formatC(bb)))

################# _________ FIXME ______________________________
ff <- "/u/maechler/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/Baharev_et_al-2017_table3.txt"
qbB2017 <- t( data.matrix(read.table(ff)) )
dimnames(qbB2017) <- dimnames(qbet)

if(inherits(p, "mpfr")) qbet <- as(qbet, "mpfr")
print(system.time(
    for(ia in seq_along(aa)) {
        a <- aa[ia]; cat("\na=",a,": b=")
        for(ib in seq_along(bb)) {
            b <- bb[ib]; cat(b," ")
            qbet[ia, ib] <- qbetaD94(p, a, b, ncp = 0, delta=delta)# def.  eps = delta^2
        }
        cat("\n")
    }
))# system.time(.)

## Safe the expensive computations !!  (FIXME: use ../inst/safe/  the same as in ~/R/Pkgs/DPQ/ !)
if(inherits(p,"mpfr")) {
    prec <- getPrec(p)
    ## FIXME  eps, delta                                    delta eps
    saveF <- paste0("qbetaD94_tab3_Rmpfr_pr", prec,
                    "_",formatC(delta),"_",formatC(eps), ".rds")
    od <- setwd("~/R/MM/Pkg-ex/Rmpfr")####################<<<<<<<<<<<<<<<<<<<<<<< FIXME >>>>>>>>>>>>>>>>>>
    saveRDS(qbet, file=saveF)
    setwd(od)
}
## number of correct digits:
summary(nDig <- asNumeric(-log10(abs(1- qbet/qbB2017))))

matplot(bb, t(asNumeric(-log10(abs(1- qbet/qbB2017)))), type = "o", log="xy")


## Looks not good:

## qbb <- qbetaD94(p, a, b, ncp = 100, delta=1e-18, verbose=3) # verbose=3: too much
## but itrmax=100'000 is not sufficient.
qbb <- qbetaD94(p, a, b, ncp = 100, delta=1e-18, itrmax = 1e8, verbose=2)
## Note how silly it is to have a very small 'eps' in a situation were 'x' is still far from the truth
## ===> Idea:  Much faster if 'eps' is "large" at the beginning, when the Newton 'd.x' will be inaccurate anyway !!

pb <- pbetaD94(p, a, b, ncp = 100)


## Other ideas:
## 1) in case Newton is not usable, be better than  x' = (1+x)/2 {on right hand} : rather use Regula Falsi, or smart unirootR() !
## 2) in these cases, use rough estimates, e.g., a few steps of unirootR()    ##    with e.g., eps=1e-3


if(FALSE) {
}
