require(DPQmpfr)
require("Rmpfr") # (is in DPQmpfr's strict dependencies)

options(warn = 1)# warnings *immediately*

##  *.time()  utilities   from
source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
source(system.file(package="DPQ", "test-tools.R", mustWork=TRUE))
## => list_() , loadList() ,  readRDS_() , save2RDS()

(doExtras <- DPQmpfr:::doExtras())
## save directory (to read from):
(sdir <- system.file("safe", package="DPQmpfr"))
## initially, when we have old version of pkg (only for pkg maintainer!!)
if(!nzchar(sdir) && dir.exists(pDir <- "~/R/Pkgs/DPQmpfr")) {
    sdir <- file.path(pDir, "inst/safe")
    if(!dir.exists(sdir)) dir.create(sdir)
}
## IGNORE_RDIFF_BEGIN
(has.sdir <- dir.exists(sdir))
if(!has.sdir)
    cat("safe directory ", sQuote(sdir), " does not exist!\n")
## IGNORE_RDIFF_END


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


qb_sfile <- function(p, precB, mpfr) {
    stopifnot(is.numeric(precB), precB == trunc(precB), precB >= 10,
              0 <= p, p < 1, is.logical(mpfr))
    paste0("tests_qbD94_", as.character(precB),"_p",
           as.character(round(100*asNumeric(p))),
           if(mpfr)"_mpfr", ".rds")
}


## Compare with  Table 3  of  Baharev_et_al 2017 %% ===> ../man/qbBaha2017.Rd <<<<<<<<<<<<
qbetaD94sim <- function(p = 0.95, # p = 1 - alpha
                        delta = NULL,
                        precBits = 128,
                        saveDir,
                        aa = c(0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 25),
                        bb = c(1:15, 10*c(2:5, 10, 25, 50)))
{
    stopifnot(dir.exists(saveDir))
    r <- lapply(c(double=FALSE, mpfr=TRUE), function(simMpfr) { ## do *both*: without / with mpfr(.)
        cat("simMpfr =", simMpfr, ":\n")
        sFile <- file.path(saveDir, qb_sfile(p, precBits, simMpfr))
        if(file.exists(sFile) && (simMpfr || !doExtras)) {
            ssR_l <- readRDS_(sFile)
            str(ssR_l)
            loadList(ssR_l, envir = environment())
            cat("[Ok]\n")
        } else { ## do run the simulation  -- always if(!simMpfr & doExtras) :
            qbet <- matrix(NA_real_, length(aa), length(bb),
                           dimnames = list(a = formatC(aa), b = formatC(bb)))
            if(!simMpfr) {
                ## p <- 0.95      # = 1 - alpha <<<--- using double precision everywhere below
                ## delta <- 1e-7  # <=> eps = delta^2 = 1e-14
            } else { ## simMpfr: takes *VERY LONG* (~ 60 minutes!) in the loop below  ((why ???))
                ## p <- 1 - 1/mpfr(20,128) # = 1 - alpha <<<--- using MPFR everywhere below
                ## delta <- 1e-18
                p <- mpfr(p, precBits=precBits)
                qbet <- as(qbet, "mpfr")
            }
            if(is.null(delta))
                delta <- if(simMpfr) 1e-18 else 1e-7
            else stopifnot(is.numeric(delta), delta >= 0)
            eps <- delta^2 # default; (delta,eps) are needed for safe-file name
            showSys.time(
                for(ia in seq_along(aa)) {
                    a <- aa[ia]; cat("\na=",a,": b=")
                    for(ib in seq_along(bb)) {
                        b <- bb[ib]; cat(b," ")
                        qbet[ia, ib] <- qbetaD94(p, a, b, ncp = 0, delta=delta)# eps := delta^2
                    } ##~~~~            ~~~~~~~~
                    cat("\n")
                }
            )
            print(summary(warnings())) # many warnings from qbeta() inaccuracies
            ## save2RDS() writes ..
            save2RDS(list_(aa, bb, qbet), file = sFile)
            cat("[Ok]\n")
        } ## --- do run simulation ---------------------------- 0
        qbet
    }) # F (simMpfr)
    ## return
    c(list_(aa, bb), r)
}

if(interactive()) ## try it out small
    rdummy <- qbetaD94sim(a=1, b=1:2, saveDir = tempdir())

if(has.sdir || doExtras)
    r.95 <- qbetaD94sim(p = 0.95, saveDir = if(has.sdir) sdir else tempdir())
str(r.95)
## $ aa
## $ bb
## $ double [qbet]
## $ mpfr   [qbet]
loadList(r.95)
rm(double, mpfr) # they shadow R functions

## number of correct digits [mpfr is only slightly better -- ??}
data(qbBaha2017, package="DPQmpfr")
dm <- with(r.95, c(length(aa), length(bb)))
stopifnot(exprs = {
    is.matrix(qbBaha2017)
    identical(dm, dim(qbBaha2017))
    identical(dm, dim(r.95$double))
    identical(dm, dim(r.95$ mpfr ))
})

ver <- setNames(, names(r.95)[3:4])
nD <- lapply(ver, function(nm) asNumeric(-log10(abs(1 - r.95[[nm]] / qbBaha2017))))

lapply(nD, round, digits=1)
## well,  Ding(1994)  may not be so good?

matplot (bb, t(nD$ double), type = "o", log="xy", xaxt="n")
matlines(bb, t(nD$ mpfr  )) # again:  mpfr is almost the same
axis(1, at=bb, cex.axis = 0.8, mgp = (2:0)/2)

showProc.time()

## Looks not good:

## NB: *relative* convergence is not good here for f() ~= 0 !!! <<< Ding *did* have absolute
## ==> use same idea as  'nls(... control = list(scaleOffset=1))' !!

a <- 1.5
b <- 7

if(interactive())##  verbose=3 now  suggest it's hopeless (?!!)
try({ ## the 'mpfr' version *does* converge  "n = 372"
    qbb <- list(mpfr = qbetaD94(mpfr(0.95, 128), shape1= a, shape2 = b, ncp = 100, delta=1e-18, itrmax = 2000, verbose=3),
                dble = qbetaD94(0.95,            shape1= a, shape2 = b, ncp = 100, delta=1e-7,  itrmax =  600, verbose=3))
    ## also itrmax = 1e6 fails the same
    ## even delta = 1e-2 fails: reason: f() = F'() = 1.25e-59 <<< bound2 ~= 3e-16
})

## Note how silly it is to have a very small 'eps' in a situation were 'x' is still far from the truth
## ===> Idea:  Much faster if 'eps' is "large" at the beginning, when the Newton 'd.x' will be inaccurate anyway !!

## Before we had 'log_scale', ncp=100  could not work
pbetaD94(.99,  a, b, ncp = 100, log_scale=FALSE, verbose=TRUE) -> pb. # now works!
## Now that pbetaD94() can use  'log_scale', this converges at n=64974 (!):
pbetaD94(.99,  a, b, ncp = 100, log_scale= TRUE, verbose=TRUE) -> pbL

all.equal(0.9999976, pb.)
all.equal(0.9999976, pbL)

showProc.time()

## Other ideas:
## 1) in case Newton is not usable, be better than  x' = (1+x)/2 {on right hand} : rather use Regula Falsi, or smart unirootR() !
## 2) in these cases, use rough estimates, e.g., a few steps of unirootR()    ##    with e.g., eps=1e-3
