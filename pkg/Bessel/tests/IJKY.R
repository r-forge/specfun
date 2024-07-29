library(Bessel)

#### Test cases for the Bessel functions I(), J(), K(), Y()
####                    -----------------------------------

all.eq <- function(x,y, tol=1e-15,...) all.equal(x, rep_len(y,length(x)), tol=tol, ...)

(op <- options(width=max(99, getOption("width")), warn = 0,
               nwarnings = 999, warnPartialMatchArgs = FALSE)) # all.equal(*,*, tol = ..)

isOSunix <- .Platform$OS.type == "unix"
windows <- !isOSunix

stop_or_w <- if(isOSunix) stop else warning

### --- For real arguments -- Comparisons with bessel[IJYK]()

x <- c(1e-6, 0.1, 1:10, 20, 100, 200)# larger x : less (relative) accuracy (??)
nus <- c(-6, -5.9, -1.1, .1, 1:4, 10, 20)
## From R 2.10.1 (or 2.10.0 patched, >= Nov.23 2009) on, this works, previously
## 	R"s  besselJ() had an inaccuracy that we have corrected here !
## FIXME ? ---- we currently need to fudge for negative nu
## note that (nu != -6 but nu ~ -6, |x| << |nu|) is still unstable,
## since sin(nu*pi) ~ 0 (but not exactly!) and besselK() is large
nS <- 10
for(F in c("I","J","K","Y")) {
    cat("\n", F," ---  nu : ")
    zF <- get(paste0("Bessel", F)) # our pkg 'Bessel'
    FF <- get(paste0("bessel", F)) # base R
    stopifnot(is.function(zF), is.function(FF))
    for(nu in nus) {
        x. <- x
        spec <- FALSE ## nu < 0 && F %in% c("I","J")
        cat(format(nu), if(spec)"* " else " ", sep="")
        zr <- zF(x., nu, nSeq = nS)
        rr <- outer(x., nu+ sign(nu)*(seq_len(nS)-1), FF)
        stopifnot(all.equal(zr, rr, tol = 500e-16))
    }; cat("\n")
}

zr0 <- if(file.exists(sf <- "zr_IJKY.rds")) { readRDS(sf)
       } else { print(saveRDS(zr, file=sf, version=2)) ; zr }

## Show:
all.equal(zr, zr0, tol = 0)

if(!isTRUE(ae <- all.equal(zr, zr0, tol = 1e-12)))
    stop_or_w(ae)


##  "limit  z -> 0  does not exist (there are many complex "Inf"s),
##		    but for z = real, z >=0 is -Inf
stopifnot(BesselY(0,1) == -Inf,# == besselY(0,1),
	  is.nan(BesselY(0+0i, 1)))

### *Large* arguments ,
### However, base::bessel*(): only I() and K() have 'expon.scaled'

x <- c(1000*(1:20), 20000*2^(1:20))
str(rI <- BesselI(x, 10, nSeq = 5, expon.scaled=TRUE))

if(getRversion() >= "2.8.2") { ## besselI(), besselJ() with larger working range
    ri2 <- outer(x, 10+seq_len(5)-1, besselI, expon.scaled=TRUE)
    stopifnot(all.equal(rI[1:20,], ri2[1:20,], tol = 8e-16))
    stopifnot(all.equal(rI[1:35,], ri2[1:35,], tol = 0.04))# base::besselI is underflowing to zero
} ## R >= 2.8.2

rI0   <- if(file.exists(sf <- "r_I.rds")) { readRDS(sf)
         } else { print(saveRDS(rI, sf, version=2)) ; rI }

ri2.0 <- if(file.exists(sf <- "r_i2.rds")) { readRDS(sf)
         } else { print(saveRDS(ri2, sf, version=2)); ri2 }

## Show the closeness (on different platforms):
all.equal(rI,  rI0,   tolerance = 0)
all.equal(ri2, ri2.0, tolerance = 0)

stopifnot(exprs = {
    all.equal(rI,  rI0,   tolerance = if(windows) 1e-10 else 1e-14)
    all.equal(ri2, ri2.0, tolerance = if(windows) 1e-10 else 1e-14)
})


## e.g. this x is too large:
x. <- 1310720000
BesselI(x., 10, nSeq = 5, expon.scaled=TRUE)
##      [,1] [,2] [,3] [,4] [,5]
## [1,]  NaN  NaN  NaN  NaN  NaN
## Warning message:
## In BesselI(1310720000, 10, nSeq = 5, expon.scaled = TRUE) :
##   'zbesi(1.31072e+09 + 0i, nu=10)'  -> ierr=4: |z| or nu too large
stopifnot(exprs = {
          (bI <- besselIasym(x., 10:14, expon.scaled=TRUE)) > 0
    all.eq(bI, bI1 <- besselIasym   (x., 10:14, expon.scaled=TRUE, k.max=1), tol=4e-15)# got 1.63e-15
    all.eq(bI, bI2 <- besselIasym   (x., 10:14, expon.scaled=TRUE, k.max=2))
    all.eq(bI, bI1n<- besselI.nuAsym(x., 10:14, expon.scaled=TRUE, k.max=1))
                                        # here, k. = 1 is slightly better
    all.eq(bI, bI2n<- besselI.nuAsym(x., 10:14, expon.scaled=TRUE, k.max=2))
    all.eq(bI, bI4n<- besselI.nuAsym(x., 10:14, expon.scaled=TRUE, k.max=4))
})

## How good are the different 'k's
library(Rmpfr)# {it has been _import_ed anyway in Bessel}
k.max <- 5
bImp <- new("mpfr", unlist(lapply(0:k.max, function(k)
    besselI.nuAsym(mpfr(x., 256), 10, expon.scaled=TRUE, k.max=k)),
    recursive = FALSE))
cbind(k = 0:(k.max-1), err.k = asNumeric(bImp - bImp[k.max+1])[-(k.max+1)])
## k      err.k
## 0 -1.051e-15
## 1 -4.510e-25
## 2 -3.584e-34
## 3 -4.187e-43
## 4 -6.469e-52


## K():  Bessel:: vs  R base:: ------------

str(rK <- BesselK(x, 10, nSeq = 5, expon.scaled=TRUE))  # our Bessel pkg [=> 20 warnings !]
rK2 <- outer(x, 10+seq_len(5)-1,
            besselK, expon.scaled=TRUE)                 # base R
stopifnot(all.equal(rK[1:35,], rK2[1:35,], tol = 8e-16))

cbind(x, local({ M <- rK; colnames(M) <- paste0("nu=",10+seq_len(ncol(M))-1); M }))


## Behaviour-------------  x --> 0  -----------------------------
## From examples in example(besselI): ~/R/D/r-devel/R/src/library/base/man/Bessel.Rd

## J():
nus <- c(0:5, 10, 20)
x0 <- 2^seq(-16, 5, length.out=256)
rbJ <- vapply(sort(c(nus, nus+0.5)), besselJ, x=x0, FUN.VALUE=x0)
rBJ <- vapply(sort(c(nus, nus+0.5)), BesselJ, z=x0, FUN.VALUE=x0)
stopifnot(all.equal(rbJ, rBJ, tol=1e-14)) # Lx 64b: 1.728e-15

## K():
x0 <- 2^seq(-10, 8, length.out=256)
rbK <- vapply(sort(c(nus, nus+0.5)), besselK, x=x0, FUN.VALUE=x0)
rBK <- vapply(sort(c(nus, nus+0.5)), BesselK, z=x0, FUN.VALUE=x0)
stopifnot(all.equal(rbK, rBK, tol=1e-14)) # Lx 64b: 1.3257e-15

## TODO?:  expon.scale  here

## Y():
x <- seq(1/32, 40, by=1/32)
rbY <- vapply(sort(c(nus, nus+0.5)), besselY, x=x, FUN.VALUE=x)
rBY <- vapply(sort(c(nus, nus+0.5)), BesselY, z=x, FUN.VALUE=x)
stopifnot(all.equal(rbY, rBY, tol=1e-14)) # Lx 64b: 1.92e-16


###--------------------- Complex  Arguments ------------------------------

besselIexpos  <- function(z, nu, expoS = TRUE) {
    drop(cbind(z,
          bI = if(is.numeric(z)) besselI(z, nu, expon.scaled=expoS) else NA
        , BI = BesselI(z, nu, expon.scaled=expoS)
        , bIa.0 = besselIasym(z, nu, k.max=0, expon.scaled=expoS)
        , bIa.1 = besselIasym(z, nu, k.max=1, expon.scaled=expoS)
        , bIa.2 = besselIasym(z, nu, k.max=2, expon.scaled=expoS)
        , bIa.3 = besselIasym(z, nu, k.max=3, expon.scaled=expoS)
        , bIa.4 = besselIasym(z, nu, k.max=4, expon.scaled=expoS)
        , bIa.6 = besselIasym(z, nu, k.max=6, expon.scaled=expoS)
        , bIa.9 = besselIasym(z, nu, k.max=9, expon.scaled=expoS)
        , bIa.19 =besselIasym(z, nu, k.max=19,expon.scaled=expoS)
        , bIna.0 = besselI.nuAsym(z, nu, k.max=0, expon.scaled=expoS)
        , bIna.1 = besselI.nuAsym(z, nu, k.max=1, expon.scaled=expoS)
        , bIna.2 = besselI.nuAsym(z, nu, k.max=2, expon.scaled=expoS)
        , bIna.3 = besselI.nuAsym(z, nu, k.max=3, expon.scaled=expoS)
        , bIna.4 = besselI.nuAsym(z, nu, k.max=4, expon.scaled=expoS)
        , bIna.5 = besselI.nuAsym(z, nu, k.max=5, expon.scaled=expoS)
          ))
}


"
z := 10000 + 10000 I
N[Exp[-Re[z]] * BesselI[10, z]]
" # -- see ../misc/MM_NUMERICS_Bessel/Bessel_I.txt
I10k1i.true <- -0.0033357343879205302021 + 0.0002661591388785316826i # from M.
bI10k1i <- besselIexpos(10000*(1+1i), nu=10) ## all look good now:
cbind(Mod(bI10k1i[-(1:2)] - I10k1i.true)/Mod(I10k1i.true))
## BI     1.619984e-17
## bIa.0  3.531183e-03
## bIa.1  6.104547e-06
## bIa.2  6.746205e-09
## bIa.3  5.233303e-12
## bIa.4  2.877969e-15
## bIa.6  2.336374e-16 < has converged:
## bIa.9  2.336374e-16
## bIa.19 2.336374e-16
## bIna.0 8.839029e-06
## bIna.1 3.525725e-10
## bIna.2 1.083583e-12
## bIna.3 1.073230e-12 -- does not get better from here ??
## bIna.4 1.073230e-12   [also not better for larger nu=200, see below]
## bIna.5 1.073230e-12

bInms <- names(bI10k1i[-(1:2)])
nIa.0 <- bInms[bInms != "bIa.0"]
n0nms <- bInms[-grep("\\.0$", bInms)]
n1nms <- n0nms[-grep("\\.1$", n0nms)]
nIhi  <- -c(1:2, 4:8, 12L)
stopifnot(exprs = {
    all.equal(    bI10k1i[["BI"]], I10k1i.true, tol=1e-15)# 1.62e-17 [Linux F28 64bit]
    all.eq(unname(bI10k1i[bInms]), I10k1i.true, tol = .0006)# 0.0002364 (L.64b)
    all.eq(unname(bI10k1i[n0nms]), I10k1i.true, tol = 1e-6)  # 4.7e-7   (L.64b)
    all.eq(unname(bI10k1i[n1nms]), I10k1i.true, tol = 2e-9)  # 6.14e-10 (L.64b)
})

bI10k1i_100 <- besselIexpos(10000*(1+1i), nu=100) ## "large nu"
bI10k1i_200 <- besselIexpos(10000*(1+1i), nu=200) ## "large nu"
## M.
I10k1i_1c.true <- -0.0025759149166597967497  - 0.0004365793917996836889i
I10k1i_2c.true <- -0.00074969703502560811294 - 0.0009802957040275765138i
## Overview ... *I.nuAsym() not really good for larger k -- "large nu" did *not* help
signif(cbind("nu=100" = Mod(bI10k1i_100[bInms] / I10k1i_1c.true - 1),
             "nu=200" = Mod(bI10k1i_200[bInms] / I10k1i_2c.true - 1)), 3)

stopifnot(exprs = {
    all.equal(    bI10k1i_100[["BI"]], I10k1i_1c.true, tol = 2e-15)# 3.55e-16 [Linux F28 64bit]
    all.eq(unname(bI10k1i_100[bInms]), I10k1i_1c.true, tol = .05  )# 0.03166  (L.64b)
    all.eq(unname(bI10k1i_100[nIa.0]), I10k1i_1c.true, tol = .01  )# 0.00596  (L.64b)
    all.eq(unname(bI10k1i_100[-(1:7)]),I10k1i_1c.true, tol = 5e-5 )# 6.55e-6  (L.64b)
    all.eq(unname(bI10k1i_100[nIhi ]), I10k1i_1c.true, tol = 8e-8 )# 1.88e-8  (L.64b)

    ## and here it's even worse [ why o why o why ??? ]
    all.equal(    bI10k1i_200[["BI"]], I10k1i_2c.true, tol = 2e-12)# 7.19e-13 [Linux F28 64bit]
    all.eq(unname(bI10k1i_200[bInms]), I10k1i_2c.true, tol = 0.5 ) # 0.325    (L.64b)  << !!!!!
    all.eq(unname(bI10k1i_200[nIa.0]), I10k1i_2c.true, tol = .01  )# 0.00596  (L.64b)
    all.eq(unname(bI10k1i_200[nIhi ]), I10k1i_2c.true, tol = .003 )# 5.99e-4  (L.64b)
})


## For now, another smaller |z| example only:
z20_5 <- 20 + 5i
I_10.z20_5 <- 0.0056200852295677786309-0.0060677028739147767132i # from M.
stopifnot(exprs = {
    all.equal(BesselI (z20_5, 10, expon.scaled=TRUE), I_10.z20_5,
              tol = 1e-15) # 2.345e-16 [Lnx_64b]
    all.equal(besselIs(z20_5, 10, expon.scaled=TRUE), I_10.z20_5,
              tol = 8e-15) # 1.049e-15 [Lnx_64b]
    ## with negative real part:
    all.equal(BesselI (-20+5i, 10, expon.scaled=TRUE), Conj(I_10.z20_5),
              tol = 1e-15) # 2.345e-16 [Lnx_64b]       ^^^^
    all.equal(besselIs(-20+5i, 10, expon.scaled=TRUE), Conj(I_10.z20_5),
              tol = 1e-14) # 2.759e-15 [Lnx_64b]       ^^^^
})
bInuAs.20_5  <- sapply(0:5, function(k.m) besselI.nuAsym(z20_5, 10, k.max=k.m))
bInuAsEX20_5 <- sapply(0:5, function(k.m) besselI.nuAsym(z20_5, 10, k.max=k.m, expon.scale=TRUE))
print(cbind(c(bInuAs.20_5, #  converging slowly to true :
              True=exp(20)*I_10.z20_5)), digits=10)
print(cbind(c(bInuAsEX20_5, #  converging slowly to true :
              True=I_10.z20_5)), digits=10)
stopifnot(exprs = {
    (err <- Mod(bInuAs.20_5*exp(-20) / I_10.z20_5 - 1)) < 6*10^-c(3,5:9)
    all.equal(err, Mod(bInuAsEX20_5 / I_10.z20_5 - 1), tol= 1e-9) # 2.5e-12
})

z0 <- round(c(c(.5, 1, 2, 5)/10, 1:10)*1000)
z <- list(1, 2-1i, 1+1i, 1-2i)
names(z) <- local({
    c <- format(lapply(z, function(.) if(Im(.)) . else Re(.)))
    ifelse(c == "1", "N", paste0("N*(",c,")")) # (no longer UTF multiplication dot)
})
z <- lapply(z, function(f) f*z0)

Iz <- lapply(z, besselIexpos, nu = 10)
## ---- for now: fine with 'numeric', not at all ok with complex !!!!
print(lapply(Iz, t), digits=4)
relE.Iz <- lapply(Iz, function(m) abs(m[,-(1:2)]/m[,"BI"] - 1))
relE.Iz ## rel.errors : now all go to "zero" nicely
## FIXME - check  Iz!


## z / 100
IzE <- lapply(lapply(z, `/`, 100), besselIexpos, nu = 10, expoS=FALSE)
## FIXME - check ... seems +- reasonable, notably the bIna.k ones !
if(FALSE)
print(lapply(IzE, t), digits=4)
relE.IzE <- lapply(IzE, function(m) abs(m[,-(1:2)]/m[,"BI"] - 1))
relE.IzE ## rel.errors nicely go to "zero" as they should
## FIXME check!

## z / 40 (<<-- barely no Inf etc)
IzE. <- lapply(lapply(z, `/`, 40), besselIexpos, nu = 10, expoS=FALSE)
if(FALSE)
print(lapply(IzE., t), digits=4)
relE.IzE. <- lapply(IzE., function(m) abs(m[,-(1:2)]/m[,"BI"] - 1))
relE.IzE. ## rel.errors nicely go to "zero" as they should
## FIXME check!


besselKexpos  <- function(z, nu, expoS = TRUE) {
    drop(cbind(z,
          bK = if(is.numeric(z)) besselK(z, nu, expon.scaled=expoS) else NA
        , BK = BesselK(z, nu, expon.scaled=expoS)
        , bKa.0 = besselKasym(z, nu, k.max=0, expon.scaled=expoS)
        , bKa.1 = besselKasym(z, nu, k.max=1, expon.scaled=expoS)
        , bKa.2 = besselKasym(z, nu, k.max=2, expon.scaled=expoS)
        , bKa.3 = besselKasym(z, nu, k.max=3, expon.scaled=expoS)
        , bKa.4 = besselKasym(z, nu, k.max=4, expon.scaled=expoS)
        , bKa.6 = besselKasym(z, nu, k.max=6, expon.scaled=expoS)
        , bKa.9 = besselKasym(z, nu, k.max=9, expon.scaled=expoS)
        , bKa.19= besselKasym(z, nu, k.max=19,expon.scaled=expoS)
        , bKna.0 = besselK.nuAsym(z, nu, k.max=0, expon.scaled=expoS)
        , bKna.1 = besselK.nuAsym(z, nu, k.max=1, expon.scaled=expoS)
        , bKna.2 = besselK.nuAsym(z, nu, k.max=2, expon.scaled=expoS)
        , bKna.3 = besselK.nuAsym(z, nu, k.max=3, expon.scaled=expoS)
        , bKna.4 = besselK.nuAsym(z, nu, k.max=4, expon.scaled=expoS)
        , bKna.5 = besselK.nuAsym(z, nu, k.max=5, expon.scaled=expoS)
          ))
}

cbind(besselKexpos(10000*(1+1i), nu=10)) # looks promising!
stopifnot(
    all.equal(
        BesselK(10000*(1+1i), nu=10, expon.scale=TRUE),
        0.0097510334110568110628-0.0040675270897257763814i, # from Mathematica
        tol=1e-15)# Linux F28 64bit: Mean relative Mod difference: 1.642e-16
)


Kz <- lapply(z, besselKexpos, nu = 10)
print(lapply(Kz, t), digits=4)
## Relative Errors [Assuming "BK" is accurate also here]: --- here all fine
relE.Kz <- lapply(Kz, function(m) abs(m[,-(1:2)]/m[,"BK"] - 1))
relE.Kz ## rel.errors nicely go to "zero" as they should
## FIXME: check these relative errors

## Open question: here the  besselK.nuAsym()  converge (to "0" ~ ke-16) for k --> 5;
#  ------------- why not for besselI.nuAsym() for the complex case [no problem in real case !!] ???


## and J() and Y()  use imaginary part for scaling by  exp( -|abs(Im(z))| ) :
stopifnot(exprs = {
    all.eq(BesselJ(20 - 5i, 7, expo=TRUE), exp(-5) * BesselJ(20 - 5i, 7))
    all.eq(BesselY(20 + 5i, 7, expo=TRUE), exp(-5) * BesselY(20 + 5i, 7))
    all.eq(BesselJ(20 + 1i, 8, expo=TRUE), exp(-1) * BesselJ(20 + 1i, 8))
    all.eq(BesselY(20 - 1i, 8, expo=TRUE), exp(-1) * BesselY(20 - 1i, 8))
})




## check Identity     J_nu(i z) =  c(nu) * I_nu(z) :

## for *integer* nu, it's simple
N <- 100
set.seed(1)
for(nu in 0:20) {
    cat(nu, "")
    z <- complex(re = rnorm(N),
                 im = rnorm(N))
    r <- BesselJ(z * 1i, nu) / BesselI(z, nu)
    stopifnot(all.equal(r, rep.int(exp(nu/2*pi*1i), N)))
}; cat("\n")

nus <- round(sort(rlnorm(20)), 2)

if(FALSE) { ## Bug ??
## For fractional nu, there's a problem (?) :
for(nu in nus) {
    cat("nu=", formatC(nu, wid=6),":")
    z <- complex(re = rnorm(N),
                 im = rnorm(N))
    r <- BesselJ(z * 1i, nu) / BesselI(z, nu)
    r.Theory <- exp(nu/2*pi*1i)
    cat("correct:"); print(table(Ok <- abs(1 - r /r.Theory) < 1e-7))
    cat("log(<wrong>) / (i*pi/2) :",
        format(unique(log(signif(r[!Ok], 6)))/(pi/2 * 1i)),
        "\n")
}

}# not for now -- Bug ?


### Replicate some testing "ideas" from  zqcbi.f (TOMS 644 test program)

##   zqcbi is a quick check routine for the complex i bessel function
##   generated by subroutine zbesi.
##
##   zqcbk generates sequences of i and k bessel functions from
##   zbesi and cbesk and checks the wronskian evaluation
##
##         I(nu,z)*K(nu+1,z) + I(nu+1,z)*K(nu,z) = 1/z
##
##   in the right half plane and a modified form
##
##        I(nu+1,z)*K(nu,zr) - I(nu,z)*K(nu+1,zr) = c/z
##
##   in the left half plane where zr: = -z and c := exp(nu*sgn*pi*i) with
##   sgn=+1 for Im(z) >= 0 and sgn=-1 for Im(z) < 0.          ^^^
##                                                          ( ||| corrected, MM)
N <- 100
nS <- 20
set.seed(257)

nus. <- unique(sort(c(nus,10*nus)))
## For exploration nus. <- (1:80)/4
for(i in seq_along(nus.)) {
    nu <- nus.[i]
    cat(nu, "")
    z <- complex(re = rnorm(N),
                 im = rnorm(N))
    P <- Re(z) >= 0
    rI  <- BesselI( z,     nu, nSeq = 1+nS) # -> for (nu, nu+1, ...,nu+nS)
    rKp <- BesselK( z[ P], nu, nSeq = 1+nS)
    rKm <- BesselK(-z[!P], nu, nSeq = 1+nS)
    ##
    sgn <- ifelse(Im(z) >= 0, +1, -1)
    Izp <- 1 / z [ P]
    for(j in 1:nS) {
        nu.. <- nu + j-1
        allEQ <- function(x,y) all.equal(x,y, tol= max(1,nu..)*nS*2e-15)
        c. <- exp(pi*nu..*sgn*1i)
        Izm <- (c./z)[!P]
        stopifnot(allEQ(rI[ P,j  ]*rKp[,j+1] + rI[ P,j+1]*rKp[,j  ], Izp),
                  allEQ(rI[!P,j+1]*rKm[,j  ] - rI[!P,j  ]*rKm[,j+1], Izm) )
    }
}; cat("\n")


### Replicate some testing "ideas" from  zqcbk.f (TOMS 644 test program)

##  part 1) in the right half plane ----> see above (I & K)
##  part 2)
##      the analytic continuation formula
##     for H(nu,2,z) in terms of the K function
##
##           K(nu,z) = c3*H(nu,2,zr) + c4*H(nu,1,zr)    Im(z) >= 0
##                   =  conjg(K(nu,conjg(z)))           Im(z) < 0
##
##     in the left half plane where c3=c1*conjg(c2)*c5, c4 = c2*c5
##     c1=2*cos(pi*nu),   c2=exp(pi*nu*i/2),   c5 =-pi*i/2   and
##     zr = z*exp(-3*pi*i/2) = z*i


### Replicate some testing "ideas" from  zqcbj.f (TOMS 644 test program)

##     zqcbj is a quick check routine for the complex J bessel function
##     generated by subroutine zbesj.
##
##     zqcbj generates sequences of J and H bessel functions from zbesj
##     and zbesh and checks the wronskians
##
##     J(nu,z)*H(nu+1,1,z) - J(nu+1,z)*H(nu,1,z) =  2/(pi*i*z)   y >= 0
##     J(nu,z)*H(nu+1,2,z) - J(nu+1,z)*H(nu,2,z) = -2/(pi*i*z)   y < 0
##
##     in their respective half planes,  y = Im(z)

N <- 100
nS <- 20
set.seed(47)

for(i in seq_along(nus.)) {
    nu <- nus.[i]
    cat(nu, "")
    z <- complex(re = rnorm(N),
                 im = rnorm(N))
    P <- Im(z) >= 0
    rJ  <- BesselJ(  z,     nu, nSeq = 1+nS) # -> for (nu, nu+1, ...,nu+nS)
    rH1 <- BesselH(1,z[ P], nu, nSeq = 1+nS)
    rH2 <- BesselH(2,z[!P], nu, nSeq = 1+nS)
    ##
    sgn <- ifelse(Im(z) >= 0, +1, -1)
    Iz <- 2/(pi*1i*z)
    for(j in 1:nS) {
        nu.. <- nu + j-1
        allEQ <- function(x,y) all.equal(x,y, tol= max(1,nu..)*nS*1e-15)
        stopifnot(allEQ(rJ[ P,j]*rH1[,j+1] - rJ[ P,j+1]*rH1[,j],  Iz[ P]),
                  allEQ(rJ[!P,j]*rH2[,j+1] - rJ[!P,j+1]*rH2[,j], -Iz[!P]) )
    }
}; cat("\n")


### TODO __FIXME__

### Replicate some testing "ideas" from  zqcby.f (TOMS 644 test program)

##     zqcby generates sequences of y bessel functions from zbesy and
##     zbesyh and compares them for a variety of values in the (z,nu)
##     space. zbesyh is an old version of zbesy which computes the y
##     function from the h functions of kinds 1 and 2.

###--------> zbesyh() in ../src/zbsubs.c  is completely unneeded (otherwise) !


##  "limit  z -> 0  does not exist (there are many complex "Inf"s),
##		    but for z = real, z >=0 is -Inf
stopifnot(BesselY(0,1) == -Inf,# == besselY(0,1),
	  is.nan(BesselY(0+0i, 1)))

## Subject: bug in Bessel package
## From: Hiroyuki Kawakatsu <...@gmail.com>, 18 Mar 2015

z <- c(0.23+1i, 1.21-1i)
nu <- -1/2
stopifnot(length(Bz.s <- BesselI(z, nu, expon.scaled=TRUE)) == length(z))
##
## Check that the exp() scaling is correct:
Bzu <- BesselI(z, nu)
stopifnot(abs(Im(sc <- Bz.s / Bzu)) < 1e-15,
	  all.eq(Re(sc), exp(-abs(Re(z)))))
## Using  nSeq > 1  -- and checking it:
options(warn=2)# warning = error {had warning of incompatible length}:
(Bz.s3 <- BesselI(z, nu, expon.scaled=TRUE, nSeq=3))
stopifnot(length(dim(Bz.s3)) == 2,
	  dim(Bz.s3) == c(length(z), 3),
	  all.eq(Bz.s3[,1], Bz.s),
	  all.eq(Bz.s3[,2], BesselI(z, nu-1, expon.scaled=TRUE)),
	  all.eq(Bz.s3[,3], BesselI(z, nu-2, expon.scaled=TRUE)))

#### From: Andrej Gajdoš <andrejg44@gmail.com>   Date: 18 Sep 2018

## originally, non-complex argument gave non-complex result for K(), Y() -- wrongly:
stopifnot(exprs = {
    all.equal(BesselK(-10, 3),
              complex(real      = -2.72527002565987e-05,
                      imaginary = -5524.11594151861),
              tol = 1e-9)
    ##	    [ tol : not sure about true accuracy]
    all.equal(BesselY(-7,2), -0.060526609468272 - 0.60283444017188i,
              tol = 1e-9)
})

### very large |x|:
## inspired from ./Gajdos-BesselK_test.R

## larger range
xL <- 10^c(0:8, 5*(2:10), 20*(3:9), 50*(4:6))
xL <- c(-rev(xL), 0, xL)
stopifnot(!is.unsorted(xL))
xs <- 1/xL # contains +Inf, fine

options(warn = 1)# show them as they occur

## base:: Bessel functions -- undefined for x < 0  as documented
bR <- cbind(x = xL,
            I = besselI(xL, 3),# many NaN, Inf + 29 warnings
            J = besselJ(xL, 3),#  ...
            K = besselK(xL, 3),
            Y = besselY(xL, 3))
bR


BesselI(-9e9, 1)# now NaN -- FIXME:  want 'Inf' (and *no* warning)!


allBessel <- function(x, nu, ...) {
    cbind(x = x,
            I = BesselI(x, nu, ...),# many NaN, Inf + 29 warnings
            J = BesselJ(x, nu, ...),#  ditto
            K = BesselK(x, nu, ...),#    "
            Y = BesselY(x, nu, ...),#    "
            H1 = BesselH(1, x, nu, ...),
            H2 = BesselH(2, x, nu, ...))
}

## This has failed with  'zbesi(-1e+300 + 0i, nu=3)' [Fortran] error ierr = 4
## then ... failed with  Error in BesselK(xL, 3)    : 'zbesk(0 + 0i, nu=3)' unexpected error 'ierr = 1'
## then ... failed with  Error in BesselH(1, xL, 3) : 'zbesh(0 + 0i, nu=3)' unexpected error 'ierr = 1'
BR3 <- allBessel(xL, 3)
options(width = 166)
BR3

(BR0 <- allBessel(xL, 0))
(BR1 <- allBessel(xL, 1))

try(## Fails:  NA/NaN/Inf in foreign function call [Inf !]
 BRs1 <- allBessel(xs, 1)
)

xs. <- sort(xs[is.finite(xs)])

(BRs0 <- allBessel(xs., 0))
(BRs1 <- allBessel(xs., 1))

(BRs3 <- allBessel(xs., 3))




options(op)
cat('Time elapsed: ', proc.time(),'\n') # for ''statistical reasons''
