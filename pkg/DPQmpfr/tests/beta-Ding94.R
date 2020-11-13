require(DPQmpfr)
require("Rmpfr") # (is in DPQmpfr's strict dependencies)

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
    pct <- summary(proc.time())# length 3, shorter names
    function(final="\n", ind=TRUE) { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- summary(proc.time())
        delta <- (pct - ot)[ind]
	##  'Time' *not* to be translated:  tools::Rdiff() skips its lines!
        cat('Time', paste0("(",paste(names(delta),collapse=" "),"):"), delta, final)
    }
})

showProc.time(,1) # 1: only 'user'

### pbeta*() : -------------------------------------------------------------
show_pbD94 <- function(...) showSys.time(show(pbetaD94(..., verbose = TRUE)))
##                                            ~~~~~~~~

show_pbD94(mpfr(0.864,  140),  5, 5, 54, eps=1e-10)# n=233 at convergence, 0.456302619295....
show_pbD94(mpfr(0.864,  140),  5, 5, 54, eps=1e-30)# n=572  "      "     , 0.45630261933697895485..

show_pbD94(mpfr(0.956,  140),  5, 5, 170, eps=1e-10)# n= 765, 0.602242164945
show_pbD94(mpfr(0.956,  140),  5, 5, 170, eps=1e-20)# n=1323, 0.60224216500116548

show_pbD94(mpfr(0.8686, 140), 10,10,  54, eps=1e-10)# n=304, 0.918779110843..
show_pbD94(mpfr(0.8686, 140), 10,10,  54, eps=1e-20)# n=497, 0.9187791109260769059699

show_pbD94(mpfr(0.8787, 140), 20,20,  54, eps=1e-10)# n=456, 0.99986765729631163
show_pbD94(mpfr(0.8787, 180), 20,20,  54, eps=1e-20)# n=

show_pbD94(mpfr(0.922,  140), 20,20, 250, eps=1e-10)# n=743, 0.964119072835964766
show_pbD94(mpfr(0.922,  180), 20,20, 250, eps=1e-20)# n=

if(doExtras) withAutoprint({
  ## MM:  Try large lambdas ==> what we need is a *RELATIVE* bound  --> now implemented
  show_pbD94(mpfr(0.93,   160), 20,25, 10000)            # n= 5175 1.008067955986..e-116
  show_pbD94(mpfr(0.93,   160), 20,25, 10000, eps=1e-20) # n= 5513 1.008067956082285281372...e-116

  show_pbD94(mpfr(0.93,   256), 20,25, 10000, eps=1e-30) # n= 5850 1.008067956082285281381936816362733357442687434229...e-116
  ## double the precBits only
  show_pbD94(mpfr(0.93,   512), 20,25, 10000, eps=1e-30) # n= 5850 1.008067956082285281381936816362733357442687434229...e-116
  ## double the precBits again + eps !!
  show_pbD94(mpfr(0.93,  1024), 20,25, 10000, eps=1e-50) # n= 6520 1.008067956082285281381936816363684457888595749378086614..e-116
  ## and much more extreme accuracy.. seems to work, too:
show_pbD94(mpfr(0.93,  1024), 20,25, 10000, eps=1e-150)# n= 9828 1.00806795608228528138193681636368445788859574937809624269896.......e-116
})

showProc.time()

### dbeta*() :

" TODO "

### qbeta*() --- careful about speed / time usage

## This now works, thanks to log_scale !
(q.5.250 <- qbetaD94(0.95, .5, 250, ncp = 0))

## but using mpfr of course works:
showProc.time()
qbetaD94(1 - 1/mpfr(20 ,64), .5, 250)                           # 0.00766110246787361539366
qbetaD94(1 - 1/mpfr(20,128), .5, 250, eps = 1e-30, delta=1e-20) # 0.0076611024678614272722717848...
showProc.time()


## Compare with  Table 3  of  Baharev_et_al 2017 %% ===> ../man/qbBaha2017.Rd <<<<<<<<<<<<
aa <- c(0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 25)
bb <- c(1:15, 10*c(2:5, 10, 25, 50))

" MM __FIXME__ :  do *both* and compare ! "
## Choose one of these
p <- 0.95      # = 1 - alpha <<<--- using double precision everywhere below
delta <- 1e-7  # <=> eps = delta^2 = 1e-14
if(FALSE) { ## This takes *VERY LONG* (> 20 minutes) in the loop below  ((why ???))
    p <- 1 - 1/mpfr(20,128) # = 1 - alpha <<<--- using MPFR everywhere below
    delta <- 1e-18
}
eps <- delta^2 # default; (delta,eps) are needed for safe-file name
qbet <- matrix(NA_real_, length(aa), length(bb),
               dimnames = list(a = formatC(aa), b = formatC(bb)))
if(inherits(p, "mpfr"))
    qbet <- as(qbet, "mpfr")
showSys.time(
    for(ia in seq_along(aa)) {
        a <- aa[ia]; cat("\na=",a,": b=")
        for(ib in seq_along(bb)) {
            b <- bb[ib]; cat(b," ")
            qbet[ia, ib] <- qbetaD94(p, a, b, ncp = 0, delta=delta)# def.  eps = delta^2
        }
        cat("\n")
    }
)

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
data(qbBaha2017, package="DPQmpfr")
summary(nDig <- asNumeric(-log10(abs(1- qbet/qbBaha2017))))

matplot(bb, t(asNumeric(-log10(abs(1- qbet/qbBaha2017)))), type = "o", log="xy", xaxt="n")
axis(1, at=bb, cex.axis = 0.8, mgp = (2:0)/2)

showProc.time()

## Looks not good:

## NB: *relative* convergence is not good here for f() ~= 0 !!! <<< Ding *did* have absolute
## ==> use same idea as  'nls(... control = list(scaleOffset=1))' !!

##  verbose=3 now  suggest it's hopeless (?!!)
try({
if(inherits(p,"mpfr")) {
    qbb <- qbetaD94(p, a, b, ncp = 100, delta=1e-18, itrmax = 2000, verbose=3)
} else {
    qbb <- qbetaD94(p, a, b, ncp = 100, delta=1e-7,  itrmax =  600, verbose=3)
    ## also itrmax = 1e6 fails the same
    ## even delta = 1e-2 fails: reason: f() = F'() = 1.25e-59 <<< bound2 ~= 3e-16
}
})

## Note how silly it is to have a very small 'eps' in a situation were 'x' is still far from the truth
## ===> Idea:  Much faster if 'eps' is "large" at the beginning, when the Newton 'd.x' will be inaccurate anyway !!

## Before we had 'log_scale', ncp=100  could not work
try( pbetaD94(.99, a, b, ncp = 100, log_scale=FALSE,
              verbose=TRUE) ) # shows 'F=NaN' ...
## Now that pbetaD94() can use  'log_scale', this converges at n=64974 (!):
pb <- pbetaD94(.99, a, b, ncp = 100, verbose=TRUE)

showProc.time()

## Other ideas:
## 1) in case Newton is not usable, be better than  x' = (1+x)/2 {on right hand} : rather use Regula Falsi, or smart unirootR() !
## 2) in these cases, use rough estimates, e.g., a few steps of unirootR()    ##    with e.g., eps=1e-3
