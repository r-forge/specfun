### <---> ../R/dgamma.R  stirlerr() definition
###       ~~~~~~~~~~~~~  ==========

## From Maple asympt(ln(GAMMA(x+1)), x, 23);
##>> ~/maple/gamma-asympt.mw
##>> ~/maple/gamma-asympt-3.tex <<

## (ln(x)-1)*x+(1/2)*ln(x)+(1/2)*ln(2*Pi)+(1/12)/x-(1/360)/x^3+(1/1260)/x^5 -
##  (1/1680)/x^7+(1/1188)/x^9-(691/360360)/x^11+(1/156)/x^13-(3617/122400)/x^15 +
## (43867/244188)/x^17-(174611/125400)/x^19+(77683/5796)/x^21 + O(1/x^23)}

require(gmp)

Sq <- c(
as.bigq(1,  12),	## 1/12
as.bigq(1,  360),	## 1/360
as.bigq(1, 1260),	## 1/1260
as.bigq(1, 1680),	## 1/1680
as.bigq(1, 1188),	## 1/1188
as.bigq(691, 360360),	## 691/360360
as.bigq(1, 156),	## 1/156
as.bigq(3617, 122400),	## 3617/122400
as.bigq(43867, 244188), ## 43867/244188
as.bigq(174611, 125400),## 174611/125400
as.bigq(77683, 5796),	## 77683/5796
NULL)

## Actually, since ~ 2021, with BernoulliQ() in {gmp}, you can directly get e.g. 16 terms :
n <- 2*(1:16)
(Sq <- abs(BernoulliQ(n))/(n*(n-1))) # == |B_n| / (n(n-1))  for n = 2,4,6,...
##
##  [1] 1/12                  1/360                 1/1260
##  [4] 1/1680                1/1188                691/360360
##  [7] 1/156                 3617/122400           43867/244188
## [10] 174611/125400         77683/5796            236364091/1506960
## [13] 657931/300            3392780147/93960      1723168255201/2492028
## [16] 7709321041217/505920

require(Rmpfr)
cbind(mpfr(Sq, 80))
## 'mpfrMatrix' of dim(.) =  (11, 1) of precision  80   bits
##  [1,]   0.083333333333333333333333368
##  [2,]  0.0027777777777777777777777776
##  [3,] 0.00079365079365079365079365075
##  [4,] 0.00059523809523809523809523806
##  [5,] 0.00084175084175084175084175104
##  [6,]  0.0019175269175269175269175262
##  [7,]  0.0064102564102564102564102561
##  [8,]   0.029550653594771241830065352
##  [9,]    0.17964437236883057316493850
## [10,]     1.3924322169059011164274315
## [11,]     13.402864044168391994478957
## [12,]     156.84828462600201730636509
## [13,]     2193.1033333333333333333318
## [14,]     36108.771253724989357173269
## [15,]     691472.26885131306710839498
## [16,]     15238221.539407416192283370


S0 <-    0.083333333333333333333333368
S1 <-   0.0027777777777777777777777776
S2 <-  0.00079365079365079365079365075 ## 1/1260
S3 <-  0.00059523809523809523809523806 ## 1/1680
S4 <-  0.00084175084175084175084175104 ## 1/1188
S5 <-   0.0019175269175269175269175262 ## 691/360360
S6 <-   0.0064102564102564102564102561 ## 1/156
S7 <-    0.029550653594771241830065352 ## 3617/122400
S8 <-     0.17964437236883057316493850 ## 43867/244188
S9 <-      1.3924322169059011164274315 ## 174611/125400
S10<-      13.402864044168391994478957 ## 77683/5796


##=======================================================================================================

## stirlerr(z)  around  z=0 :
## ---------------------------  this is really a very different series:
# In Wolfram Alpha  :
'
Series[LogGamma[z+1] - (z+1/2)*Log[z] + z - Log[2 Pi]/2, {z, 0, 12}]
'
## ==>
"https://www.wolframalpha.com/input?i=Series%5BLogGamma%5Bz%2B1%5D+-+%28z%2B1%2F2%29*Log%5Bz%5D+%2B+z+-+Log%5B2+Pi%5D%2F2%2C+%7Bz%2C+0%2C+5%7D%9D"

## output, copyable text:
-1/2 log(2 π z) + z (-log(z) - gamma + 1) + (π^2 z^2)/12 + (z^3 polygamma(2, 1))/6 + (π^4 z^4)/360 + (z^5 polygamma(4, 1))/120 + (π^6 z^6)/5670 + (z^7 polygamma(6, 1))/5040 + (π^8 z^8)/75600 + (z^9 polygamma(8, 1))/362880 + (π^10 z^10)/935550 + (z^11 polygamma(10, 1))/39916800 + (691 π^12 z^12)/7662154500 + O(z^13)
## (generalized Puiseux series)


## MM: R code
-1/2 log(2 π z) + z (-log(z) - gamma + 1) + (π^2 z^2)/12 + (z^3 psigamma(1, 2))/6 + (π^4 z^4)/360 + (z^5 psigamma(1, 4))/120 + (π^6 z^6)/5670 + (z^7 psigamma(1, 6))/5040 + (π^8 z^8)/75600 + (z^9 psigamma(1, 8))/362880 + (π^10 z^10)/935550 + (z^11 polygamma(10, 1))/39916800 + (691 π^12 z^12)/7662154500 ## + O(z^13)
## (generalized Puiseux series)

setwd('/u/maechler/R/D/R-forge/specfun/pkg/DPQ/Misc/')
system("eog stirlerr_z=0_series_alpha.png &")

## NOTA BENE:  polygamma(der, z) == psigamma(z, deriv)  is
## ---------:  *and* polygamma(2k, 1) =  psigamma(1, 2k) ==  - (2k)! zeta(2k+1)
##                                       ~~~~~~~~~~~~~~~     ~~~~~~~~~~~~~~~~~~
## *and* fortunately,  zeta() exists in Rmpfr !!

## z = 0 above has coefficients
## c_{2k} = \psigamma(1, 2k) / (2k+1)!
          ## = - (2k)! zeta(2k+1) / (2k+1)! =
          ## = - zeta(2k+1) / (2k+1)


##--- the same, series around z = 5: ---------------------------------------
##                            ~~~~~~
'
Series[LogGamma[z+1] - (z+1/2)*Log[z] + z - Log[2 Pi]/2, {z, 5, 9}]
'
system("eog stirlerr_z=5_series_alpha.png &")

## output, copyable text:
(5 - (11 log(5))/2 + log(120) - 1/2 log(2 π)) + (z - 5) (131/60 - gamma - log(5)) + (π^2/12 - 5917/7200) (z - 5)^2 + (z - 5)^3 (2/375 + polygamma(2, 6)/6) + (π^4/360 - 14025553/51840000) (z - 5)^4 + (z - 5)^5 (3/62500 + polygamma(4, 6)/120) + (π^6/5670 - 47464869601/279936000000) (z - 5)^6 + (z - 5)^7 (1/1640625 + polygamma(6, 6)/5040) + (π^8/75600 - 1180524695078791/9405849600000000) (z - 5)^8 + (z - 5)^9 (1/140625000 + polygamma(8, 6)/362880) + O((z - 5)^10)
## (Taylor series)

## output, "Wolfram Language" (Mathematica)  text:
(5 - (11 Log[5])/2 + Log[120] - Log[2 Pi]/2) + (131/60 - EulerGamma - Log[5]) (-5 + z) + (-5917/7200 + Pi^2/12) (-5 + z)^2 + (2/375 + PolyGamma[2, 6]/6) (-5 + z)^3 + (-14025553/51840000 + Pi^4/360) (-5 + z)^4 + (3/62500 + PolyGamma[4, 6]/120) (-5 + z)^5 + (-47464869601/279936000000 + Pi^6/5670) (-5 + z)^6 + (1/1640625 + PolyGamma[6, 6]/5040) (-5 + z)^7 + (-1180524695078791/9405849600000000 + Pi^8/75600) (-5 + z)^8 + (1/140625000 + PolyGamma[8, 6]/362880) (-5 + z)^9 + O[-5 + z]^10

## MM: R-code s/ polygamma(k,z) / psigamma(z, k) /
(5 - (11 log(5))/2 + log(120) - 1/2 log(2 π)) + (z - 5) (131/60 - gamma - log(5)) + (π^2/12 - 5917/7200) (z - 5)^2 + (z - 5)^3 (2/375 + psigamma(6, 2)/6) + (π^4/360 - 14025553/51840000) (z - 5)^4 + (z - 5)^5 (3/62500 + psigamma(6, 4)/120) + (π^6/5670 - 47464869601/279936000000) (z - 5)^6 + (z - 5)^7 (1/1640625 + psigamma(6, 6)/5040) + (π^8/75600 - 1180524695078791/9405849600000000) (z - 5)^8 + (z - 5)^9 (1/140625000 + psigamma(6, 8)/362880) ##+ O((z - 5)^10)

require(Rmpfr)

## In all cases,  psigamma(2k, z+1) = rational  - (2k)! zeta(2k+1)
##                ~~~~~~~~~~~~~~~~~   ........     ~~~~~~~~~~~~~~~~   where 'rational = 0' for z = 0
kk0 <- 4L
twok1 <- mpfr(twok1, 128)
str(twok1 <- 2L*(1:kk0) + 1L) # 2k + 1
##        - (2k)!          zeta(2k+1)
F2zeta <- - gamma(twok1) * zeta(twok1)
## \psi^{(2k)}(1) =
(psigam_1 <- F2zeta)
## [1] -2.404113806319188570799476323022899981535 -24.88626612344087823195277167496882003344
## [3] -726.0114797149844353246542358918536669099  -40400.9783987476348853278236554508542788

## psigam_1 <- -(twok1-1)! zeta(twok1)
(psigam_1 <- - gamma(twok1) * zeta(twok1))

## output, copyable text:
(4 - 9 log(2) + log(24) - 1/2 log(2 π)) + (z - 4) (47/24 - gamma - log(4)) + 1/576 (48 π^2 - 473) (z - 4)^2 + (z - 4)^3 (1/128 + polygamma(2, 5)/6) + (π^4/360 - 44873/165888) (z - 4)^4 + ((z - 4)^5 (3 + 256 polygamma(4, 5)))/30720 + (π^6/5670 - 30376837/179159040) (z - 4)^6 + (z - 4)^7 (1/688128 + polygamma(6, 5)/5040) + (π^8/75600 - 6044274287/48157949952) (z - 4)^8 + ((z - 4)^9 polygamma(8, 5))/362880 + O((z - 4)^10)
##(Taylor series)
## MM: R-code s/ polygamma(k,z) / psigamma(z, k) /
(4 - 9 log(2) + log(24) - 1/2 log(2 π)) + (z - 4) (47/24 - gamma - log(4)) + 1/576 (48 π^2 - 473) (z - 4)^2 + (z - 4)^3 (1/128 + psigamma(5, 2)/6) + (π^4/360 - 44873/165888) (z - 4)^4 + ((z - 4)^5 (3 + 256 psigamma(5, 4)))/30720 + (π^6/5670 - 30376837/179159040) (z - 4)^6 + (z - 4)^7 (1/688128 + psigamma(5, 6)/5040) + (π^8/75600 - 6044274287/48157949952) (z - 4)^8 + ((z - 4)^9 psigamma(5, 8))/362880 ## + O((z - 4)^10)


## z = 4
## ===== --------------------------------
'                                                            #
Series[LogGamma[z+1] - (z+1/2)*Log[z] + z - Log[2 Pi]/2, {z, 4, 9}]
'
system("eog stirlerr_z=4_series_alpha.png &")
## output, copyable text:
(4 - 9 log(2) + log(24) - 1/2 log(2 π)) + (z - 4) (47/24 - gamma - log(4)) + 1/576 (48 π^2 - 473) (z - 4)^2 + (z - 4)^3 (1/128 + polygamma(2, 5)/6) + (π^4/360 - 44873/165888) (z - 4)^4 + ((z - 4)^5 (3 + 256 polygamma(4, 5)))/30720 + (π^6/5670 - 30376837/179159040) (z - 4)^6 + (z - 4)^7 (1/688128 + polygamma(6, 5)/5040) + (π^8/75600 - 6044274287/48157949952) (z - 4)^8 + ((z - 4)^9 polygamma(8, 5))/362880 + O((z - 4)^10)
## (Taylor series)
## MM: R-code s/ polygamma(k,z) / psigamma(z, k) /
(4 - 9 log(2) + log(24) - 1/2 log(2 π)) + (z - 4) (47/24 - gamma - log(4)) + 1/576 (48 π^2 - 473) (z - 4)^2 + (z - 4)^3 (1/128 + psigamma(5, 2)/6) + (π^4/360 - 44873/165888) (z - 4)^4 + ((z - 4)^5 (3 + 256 psigamma(5, 4)))/30720 + (π^6/5670 - 30376837/179159040) (z - 4)^6 + (z - 4)^7 (1/688128 + psigamma(5, 6)/5040) + (π^8/75600 - 6044274287/48157949952) (z - 4)^8 + ((z - 4)^9 psigamma(5, 8))/362880 ## + O((z - 4)^10)


## z = 4  ==> evaluate psigamma^(2k)(z+1) = psigamma(z+1, 2*k)
## ===== --------------------------------
## Polygamma(2,5)
##    FunctionExpand[PolyGamma[2, 5]]
'2 (2035/1728 - ζ(3))'
2*(BZ(2035)/1728 - Rmpfr::zeta(3))
psigamma(5, deriv=2) # -0.04878973

## Polygamma(4,5)
##    FunctionExpand[PolyGamma[4, 5]]
'24 (257875/248832 - ζ(5))'
psigamma(5, deriv=4) # -0.01406319

## Polygamma(6,5)
##    FunctionExpand[PolyGamma[6, 5]]
'720 (36130315/35831808 - ζ(7))' # -0.013316295488550550991684451524448804905878833309317180114640224
psigamma(5, deriv=6) # -0.0133163

## Polygamma(8,5)
##    FunctionExpand[PolyGamma[8, 5]]
'40320 (5170139875/5159780352 - ζ(9))' #-0.026121932577157847256668711028032780542422983569899521811529705...
psigamma(5, deriv=8) #  -0.02612193

##===========================================================================
##----------------- R code -------------------------
Const("gamma", 128) #  'mpfr' number ..  128 bits -- 0.5772156649015328606065120900824024310432
Const("pi",    128)

## no Rmpfr here, for now:
gamma <- 0.5772156649015328606065120900824024310432
π     <- 3.141592653589793238462643383279502884195
## z = 0:
f0. <- function(z)
-1/2*log(2*π*z) + z*(-log(z) - gamma + 1) + (π^2 * z^2)/12 + (z^3* psigamma(1, 2))/6 + (π^4 * z^4)/360 + (z^5* psigamma(1, 4))/120 + (π^6 * z^6)/5670 + (z^7* psigamma(1, 6))/5040 + (π^8 * z^8)/75600 + (z^9* psigamma(1, 8))/362880 + (π^10 * z^10)/935550 + (z^11 * psigamma(1, 10))/39916800 + (691 * π^12 * z^12)/7662154500 ## + O(z^13)
##==> see f0() further below, after f00()

## z = 4:
f4. <- function(z)
(4 - 9*log(2) + log(24) - 1/2*log(2*π)) + (z - 4) * (47/24 - gamma - log(4)) + 1/576*(48*π^2 - 473) * (z - 4)^2 + (z - 4)^3*(1/128 + psigamma(5, 2)/6) + (π^4/360 - 44873/165888) * (z - 4)^4 + ((z - 4)^5*(3 + 256* psigamma(5, 4)))/30720 + (π^6/5670 - 30376837/179159040) * (z - 4)^6 + (z - 4)^7*(1/688128 + psigamma(5, 6)/5040) + (π^8/75600 - 6044274287/48157949952) * (z - 4)^8 + ((z - 4)^9* psigamma(5, 8))/362880 ## + O((z - 4)^10)
## improve:
Pi <- Const("pi", 128)
c0 <- log(2*Pi)/2
L2 <- log(mpfr(2, 128))
L3 <- log(mpfr(3, 128))
Gamma <- Const("gamma", 128)
gamma <- asNumeric(Gamma)
BZ <- gmp::as.bigz

## really too crude approximation
f000 <- function(z) -log(z)/2
x <- 10^-(308:80)
relE000 <- asNumeric(relErrV(stirlerrM(mpfr(x, 256)), f000(x)))
summary(relE000)
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
## 0.002598 0.003190 0.004131 0.004752 0.005860 0.010078
##
#  ===> far away from accurate ... forget about it !

## This is good however for small z :
C0 <- asNumeric(c0) # log(2*Pi)/2 ~= 0.9189
f00 <- function(z) -(log(z)/2 + C0)

str(x <- lseq(2^-100, 2^-10, length=1000))
plot(x, asNumeric(relErrV(stirlerrM(mpfr(x, 256)), f00(x))), log="x", type="l")
abline(h=c(-2:2)*2^-53, lty=3, col="gray")
## rel.err <= 0.0025 : "usable" for plotting  but not as approx
##
## ==> much smaller 'x' :

x <- 10^-seq(30,14, length=1000)
x <- 10^-seq(18,15.5, length=1000)
relE00 <- asNumeric(relErrV(stirlerrM(mpfr(x, 256)), f00(x)))

plot(x,   relE00, type="l", log="x", ylim = c(-1,1)*1e-15,
     main = "relative Error of  -(log(x)/2 + C0)  approx. to  stirlerr(x)")
abline(h=c(-2:2)*2^-53, lty=3, col="gray")

plot(x,  abs(relE00), type="l", log="xy",
     main = "|relative Error|  of  -(log(x)/2 + C0)  approx. to  stirlerr(x)")
lines(x, abs(relE000), col = adjustcolor(4, 1/2), lwd=2)
lines(lowess(x, abs(relE00), f=.1), col=2, lwd=2)
abline(h=c(1:2,4)*2^-53, lty=3, col="gray")
p2 <- -(57:54)
abline(v=2^p2, col="lightblue", lty=2)
axis(1, at=2^p2, labels=as.expression(lapply(p2+0, \(n) substitute(2^N, list(N=n)))), lwd=NA, line=-2)
## ==>  f00() is  good for n < 1e-17 or  n <= 2^{-56} = 1.3878e-17
##      ~~~~~                            ~~~~~~~~~~~~

kk0 <- 4L
str(twok1 <- 2L*(1:kk0) + 1L) # 2k + 1
twok1 <- mpfr(twok1, 128)
##        - (2k)!          zeta(2k+1)
## F2zeta <- - gamma(twok1) * zeta(twok1)
## ## \psi^{(2k)}(1) =
## (psigam_1 <- F2zeta)
## [1] -2.404113806319188570799476323022899981535 -24.88626612344087823195277167496882003344
## [3] -726.0114797149844353246542358918536669099  -40400.9783987476348853278236554508542788

## psigam_1 <- -(twok1-1)! zeta(twok1)
psigam_1 <- - gamma(twok1) * zeta(twok1)

c01  <- asNumeric(1 - Gamma)
c02  <- asNumeric(Pi^2 / 12)
c0z1 <- asNumeric(-zeta(3)/3) # == psigamma(1, 2)/6
c03  <- asNumeric(Pi^4 / 360)
c0z2 <- asNumeric(-zeta(5)/5) # == psigamma(1, 4)/120
c04  <- asNumeric(Pi^6 / 5670)
c0z3 <- asNumeric(-zeta(7)/7) # == psigamma(1, 6)/5040
c05  <- asNumeric(Pi^8/75600)
c0z4 <- asNumeric(-zeta(9)/9) # == psigamma(1, 8)/362880
c06  <- asNumeric(Pi^10/935550)
c0z5 <- asNumeric(-zeta(11)/11) # == psigamma(1, 10)/39916800
c07  <- asNumeric(691 * Pi^12/7662154500)

## f0. <- function(z)
## -1/2*log(2*π*z) + z*(-log(z) - gamma + 1) + (π^2 * z^2)/12 + (z^3* psigamma(1, 2))/6 + (π^4 * z^4)/360 + (z^5* psigamma(1, 4))/120 + (π^6 * z^6)/5670 + (z^7* psigamma(1, 6))/5040 + (π^8 * z^8)/75600 + (z^9* psigamma(1, 8))/362880 ## + O(z^10)

## FYI   f00 <- function(z) -(log(z)/2 + C0)

f0 <- function(z) {
    lz <- log(z)
    -(lz/2 + C0) + # = -(log(z) + log(2 pi))/2 = - log(2 pi z)/2
        z*(-lz + c01 +
           z*(c02 + # π^2/12
              z*(c0z1 + # psigamma(1, 2)/6
                 z*(c03 + # π^4/360 +
                    z*(c0z2 + # psigamma(1, 4)/120
                       z*(c04 + # π^6/5670
                          z*(c0z3 + # psigamma(1, 6)/5040
                             z*(c05 + # π^8/75600
                                 z*(c0z4 + # psigamma(1, 8)/362880
                                    + z*(c06 + # π^10/935550
                                         z*(c0z5 + # psigamma(1, 10)/39916800
                                            z*c07 # (691 π^12)/7662154500   + O(z^13)
                                         )))))))))))
}

f0.ord.max <- 12L # and below
f0_ <- function(z, order) {
    stopifnot(0L <= (order <- as.integer(order)), order <= f0.ord.max,
              length(order) == 1L)
    lz <- log(z)
    f00 <- -(lz/2 + C0)
    switch(order + 1L, # so we start with 1 {for "integer" switch()}
           f00, # 0
           f00 + z*(-lz + c01), # 1
           f00 + z*(-lz + c01 + z* c02), # 2
           f00 + z*(-lz + c01 + z*(c02 + z* c0z1)), # 3
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z* c03))), # 4
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z* c0z2)))), # 5
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z*(c0z2 + z* c04))))), # 6
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z*(c0z2 + z*(c04 + z* c0z3)))))), # 7
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z*(c0z2 + z*(c04 + z*(c0z3 + z* c05))))))), # 8
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z*(c0z2 + z*(c04 + z*(c0z3 + z*(c05 + z* c0z4)))))))), # 9
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z*(c0z2 + z*(c04 + z*(c0z3 + z*(c05 + z*(c0z4 + z* c06))))))))), # 10
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z*(c0z2 + z*(c04 + z*(c0z3 + z*(c05 + z*(c0z4 + z*(c06 + z* c0z5)))))))))), # 11
           f00 + z*(-lz + c01 + z*(c02 + z*(c0z1 + z*(c03 + z*(c0z2 + z*(c04 + z*(c0z3 + z*(c05 + z*(c0z4 + z*(c06 + z*(c0z5 + z* c07))))))))))) # 12
           )
}


1


## (4 - 9*log(2) + log(24) - 1/2*log(2*π)) - accurately
c4.0 <- asNumeric(4 - 6*L2 + L3 - c0)
c4.1 <- asNumeric(BZ(47)/24 - gamma - 2*L2) # 47/24 - gamma - log(4)
c4.2 <- asNumeric((48*Pi^2 - 473)/576) # (48*π^2 - 473)/576  = 0.00128647..
c4.3 <- asNumeric(1/128 + (BZ(2035)/1728 - Rmpfr::zeta(3))/3) # (1/128 + psigamma(5, 2)/6)
c4.4 <- asNumeric(Pi^4/360 - BZ(44873)/165888) # π^4/360 - 44873/165888
c4.5 <- asNumeric((3 + 256 * 24*(BZ(257875)/248832 - Rmpfr::zeta(5)))/30720) # (3 + 256* psigamma(5, 4))/30720
c4.6 <- asNumeric(Pi^6/5670 - BZ(30376837)/179159040) # π^6/5670 - 30376837/179159040

f4.2 <- function(z) {
    z_4 <- z - 4
    c4.0 +
        z_4 *(c4.1 + # (47/24 - gamma - log(4))
              z_4 *((48*π^2 - 473)/576 +
                    (1/128 + psigamma(5, 2)/6)*z_4 +
                    (π^4/360 - 44873/165888) * z_4^2 +
                    (3 + 256* psigamma(5, 4))/30720 * z_4^3 +
                    (π^6/5670 - 30376837/179159040) * z_4^4 +
                    (1/688128 + psigamma(5, 6)/5040)* z_4^5 +
                    (π^8/75600 - 6044274287/48157949952) * z_4^6 +
                    psigamma(5, 8)/362880 * z_4^7)) ## + O((z - 4)^10)
}

f4 <- function(z) {
    z_4 <- z - 4
    c4.0 +
        z_4 *(c4.1 +
              z_4 *(c4.2 +
                    z_4 *(c4.3 +
                          z_4 *(c4.4 +
                                z_4 *(c4.5 +
                                      z_4 *(c4.6 +
                                            z_4 *((1/688128 + psigamma(5, 6)/5040) +
                                                  z_4 *((π^8/75600 - 6044274287/48157949952) +
                                                        z_4 *(psigamma(5, 8)/362880))))))))) ## + O((z - 4)^10)
}


## z = 5:
f5 <- function(z)
(5 - (11*log(5))/2 + log(120) - 1/2*log(2*π)) + (z - 5) * (131/60 - gamma - log(5)) + (π^2/12 - 5917/7200) * (z - 5)^2 + (z - 5)^3*(2/375 + psigamma(6, 2)/6) + (π^4/360 - 14025553/51840000) * (z - 5)^4 + (z - 5)^5*(3/62500 + psigamma(6, 4)/120) + (π^6/5670 - 47464869601/279936000000) * (z - 5)^6 + (z - 5)^7*(1/1640625 + psigamma(6, 6)/5040) + (π^8/75600 - 1180524695078791/9405849600000000) * (z - 5)^8 + (z - 5)^9*(1/140625000 + psigamma(6, 8)/362880) ##+ O((z - 5)^10)


## try/test  them ===============================================================

## 1.  z = 0 :
f0xy <- curve(f0, 0, 1.1, n = 1001)
lines(stirlerr(x) ~ x, data=f0xy,  lwd = 3, col = adjustcolor(2, 1/2))
                                        # "looks good" ...
                                        # but how large is the error?
## even only in [0, 1/4] not good enough:
## x <- seq(0, 0.15, length=2000)
x <- lseq(1e-20, 0.15, length=2000)
## more zoom in:
x <- lseq(1e-8,  0.11, length=2000)
xM <- mpfr(x, 256)
relE0  <- asNumeric(f0 (x) / stirlerrM(xM) -1)
relE0. <- asNumeric(f0.(x) / stirlerrM(xM) -1)

plot (x, abs(relE0), type = "l", lwd=2, log="xy", axes=FALSE, frame=TRUE)# fine
eaxis(1, nintLog=20); eaxis(2, nintLog = 18)
lines(x, abs(relE0.), col=adjustcolor(2, 1/2), lwd=2) # f0.() is a bit worse than f0()  <==> Horner-like nesting *does* pay
## f0.() is slightly worse in the "interesting" region [10^-8, 10^-1
abline(h=c(1,2,4)*2^-53, lty=3, col="gray")
p2 <- -(9:3)
abline(v=2^p2, col="lightblue", lty=2)
axis(1, at=2^p2, labels=as.expression(lapply(p2+0, \(n) substitute(2^N, list(N=n)))), lwd=NA, line=-2)
## ===> f0() covers for all z = n <= 0.03 (w/ O(z^{10})) and 0.07 (w/ O(z^{13})
##      ~~~~                ~~~~~~~~~~~~~                    ====

###==== all different orders:  ============= f0_(x, order) ===============================================
x <- lseq(1e-20,  0.11, length=2500)
xM <- mpfr(x, 512)
system.time( ## takes time
relEall <- vapply(0:f0.ord.max,
                  function(io) asNumeric(f0_(x, io) / stirlerrM(xM) - 1), x)
) # ~ 6.55 sec (7 sec on lynne) elapsed

matplot(x, relEall, type="l", log="x", ylim = c(-1,1)*1e-13, axes=FALSE)
eaxis(1, nintLog=20); eaxis(2) #                      ^^^^^    ==> zoom in 100 x :

## 13:1 : plot in reverse order
cols <- adjustcolor(13:1, 3/4)
ltys <- 13:1

## rel.error
matplot(x, relEall[,13:1], type="l", log="x", ylim = c(-1,1)*1e-15, axes=FALSE,
        col=cols, lty=ltys,
        xlab=quote(n), main = "stirlerr(n) series approx at n=0, order 0:12")
eaxis(1, nintLog=20); eaxis(2)
abline(h=c(-2:2)*2^-53, lty=c(3,3,1,3,3), lwd=c(1,1,2,1,1),
       col=adjustcolor("gray20", 1/2))
legend("topleft", legend=paste("order", 0:12, sep=" = "), ncol = 5,
       col=rev(cols), lty=rev(ltys), bty="n")


## |rel.error|
matplot(x, abs(relEall[,13:1]), type="l", log="xy", ylim = c(1, 100)*1e-17, axes=FALSE,
        xlab=quote(n), col=cols, lty=ltys,  ylab = "| rel.Err |",
        main = "stirlerr(n) series approx at n=0, order 0:12")
eaxis(1, nintLog=20); eaxis(2, cex.axis = .8)
abline(h = c(1,2,4)*2^-53, lty=3, col=adjustcolor("gray20", 1/2))
legend("topleft", legend=paste("order", 0:12, sep=" = "), ncol = 5,
       col=rev(cols), lty=rev(ltys), bty="n")

##====================================================

## 2.  z = 4 :
## -----------
f4xy <- curve(f4, 3, 5, n = 1001)
lines(stirlerr(x) ~ x, data=f4xy,  lwd = 3, col = adjustcolor(2, 1/2))
                                        # looks good ... but how large is the error :

##
x <- seq(3.5, 4.5, by = 2^-10)
xM <- mpfr(x, 256)
plot (x, abs(asNumeric(f4(x)/stirlerrM(xM) -1)), type = "l", lwd=2, log="y")
lines(x, abs(asNumeric(f4.2(x)/stirlerrM(xM) -1)), col=2)
lines(x, abs(asNumeric(f4. (x)/stirlerrM(xM) -1)), col=3)
abline(h = 2^-53*1:4, lty=3)
## good enough (53-bit) in ca [3.9, 4.1].....

################## Use Chebyshev - Approximation ###############


## 1. -- the full relevant interval, say [1/16, 8]

require(Rmpfr); require(DPQmpfr) # for stirlerrM()
require(sfsmisc)

x <- lseq(1/16, 8, length=2000)
xM <- mpfr(x, 256)
stirle <- stirlerrM(xM)
Pi <- Const("pi", 256)
C0 <- log(2*Pi)/2

## nicer function: "subtract" the main term (close x=0) of stirlerr
str(stirle1 <- stirle + log(2 * Pi * xM)/2)
## Class 'mpfr' [package "Rmpfr"] of length 2000 and precision 256
##  0.202829694846 0.203177914207 ...

plot(x, asNumeric(stirle1), type="l", log="x")
## does look "simple" nice (in  log(x) - scale !!)

require(pracma)

stirle1 <- function(x, prec = 256) {
    xM <- mpfr(x, prec)
    asNumeric(stirlerrM(xM) + log(xM)/2)
}

## 1 "half" term more:
stirle2 <- function(x, prec = 256) {
    xM <- mpfr(x, prec)
    asNumeric(stirlerrM(xM) + log(xM)*(1/2 + xM))
}

str(s1L <- curve(stirle1, 1/16, 8, log="x", n=1000))
##                     ^
str(s2L <- curve(stirle2, 1/16, 8, log="x", n=1000)) # hmm, not better

plot(y ~ x, data=s1L, type = "l")
## the one I want to see
plot(y ~ x, data=s1L, type = "l", log="x")
## the one to base the interpolation on --- same curve; different x-axis labels/scale
plot(y ~ log(x), data=s1L, type = "l")

## stirle2 :
plot(y ~ x, data = s2L, type = "l", log="x")

require(pracma)

chebApprox # defines the polynomial cP and coefficients for fun;
           # then evaluates these at x:
function (x, fun, a, b, n)
{
    cP <- chebPoly(n)
    cC <- chebCoeff(fun, a, b, n)
    p <- drop(cC %*% cP)
    c0 <- cC[1]
    xx <- (2 * x - (b + a))/(b - a)
    yy <- polyval(p, xx) - c0/2
    return(yy)
}

## on log2-scale  {x, a, b} ;
## 2^a = 2^-4 = 1/16
## 2^b = 2^3  =  8
l2x <- seq(-4, 3, length=1000) # log2(x)

ya50 <- chebApprox(x = l2x, fun = \(lx) stirle1(2^lx),
                   a= -4, b = 3,
                   n = 50)

str(s1L <- curve(stirle1, 1/16, 8, log="x", n=200, lwd = 3, col=2))
lines(2^l2x, ya50)

str(relE50 <- relErrV(stirle1(2^l2x), ya50))
summary(relE50)
plot(2^l2x,     relE50,  log="x", type="l")

plot(2^l2x, abs(relE50), log="xy", type="l")
abline(h = 2^-53*1:4, lty=3) # too large

## 100 chebychev coefficients {500: really all "beyond" -- because it did *not* to 'TODO"/
ya100 <- chebApprox(x = l2x, fun = \(lx) stirle1(2^lx),
                   a= -4, b = 3,
                   n = 100)
## completely diverges at boundaries !  --->  pracma's TODO ???

str(relE100 <- relErrV(stirle1(2^l2x), ya100))
summary(relE100)
plot(2^l2x,     relE100,  log="x", type="l")

plot(2^l2x, abs(relE100), log="xy", type="l")
abline(h = 2^-53*1:4, lty=3) # too large

## NB: our DPQ  chebyshevEval() *does* use Rmathlib's C code *with* the Clenshaw algo
## NB2: we *should be able to write  chebApprox()
## ---  learning from pracma's  chebCoeff() and using our DPQ
