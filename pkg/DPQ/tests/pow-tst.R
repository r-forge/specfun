require("DPQ")
require("Rmpfr")

(doExtras <- DPQ:::doExtras())
(noLdbl <- (.Machine$sizeof.longdouble <= 8)) ## TRUE when --disable-long-double

k1 <- 44300 + (1:800)
p1 <- (63/64)
x <- p1; stopifnot(  (1-x)-1 == -x , 1-x == 2^-6 )
px <- p1 ^ k1
dput(         px,, "hex")
dput(2^1000 * px,, "hex")


x <- p2 <- 1023/1024
stopifnot(  (1-x)-1 == -x , 1-x == 2^-10 )

k2 <- 700000 + (1:1000)
dput(p2^ k2, control = "digits17") ## last 2-3 digits printed *differ*
dput(p2^ k2, control = "hex")

x <- p3 <- (1- 2^-14)
stopifnot(  (1-x)-1 == -x , 1-x == 2^-14 )
log(1e-300)/ log(p3) # 11317321

k3 <- 11321012 + (1:1000)
dput(p3^ k3, control = "digits17") ## last 2-3 digits printed *differ*
dput(p3^ k3, control = "hex")

require(sfsmisc) # for relErrV()
## relative error of  x ^ k  computation - comparing with truth = <x_mpfr_{precBits])
relEP <- function(x, k, precBits) { # see below about  precBits needed
    stopifnot(is.numeric(x))
    asNumeric(relErrV(mpfr(x, precBits)^k, x^k))
}

stopifnot(.Machine$double.xmin == 2^-1022) # will use  2^-1022
osV <- abbreviate(sub("\\(.*", "", osVersion), 10)

re1. <- relEP(p1, k1, precBits = 2^12) # high precision
re1  <- relEP(p1, k1, precBits = 256)
stopifnot(re1. == re1) # ===> 256 bit suffice here
summary(re1)
summary(abs(re1))

if(!dev.interactive(TRUE)) pdf(paste0("pow-ex_", osV, ".pdf"),
                               width = 9, height=5)
plot(k1, re1, type="l", main = "rel.error {wrt MPFR} of  (63/64)^k")
mtext(sfsmisc::shortRversion(spaces=FALSE))
abline(v = min(k1[p1^k1 < 2^-1022]), lty=3)
ll <- c(1,1,2,1,1)
abline(h=(-2:2)*2^-53, lty=3-ll, lwd=ll, col=adjustcolor(ll, 1/2))
## looks perfect (in Linux; *not* windows

re2  <- relEP(p2, k2, 2^ 9)
re2. <- relEP(p2, k2, 2^14)
stopifnot(re2. == re2)
summary(re2)
summary(abs(re2))
str(rat2 <- c(MASS:::.rat(p2) $ rat))

plot(k2, re2, type="l", xlab = quote(k), main = substitute("rel.Err of " ~ (N/Z)^k,
                                                           list(N=rat2[1], Z=rat2[2])),
     ylim = range(0, re2))
if(any((px <- p2^k2) < 2^-1022)) abline(v = min(k2[px < 2^-1022]), lty=3)
abline(h=(-2:2)*2^-53, lty=3-ll, lwd=ll, col=adjustcolor(ll, 1/2))

re3 <- relEP(p3, k3, prec = 256)
summary(re3)
str(rat3 <- c(MASS:::.rat(p3, max.denominator=1e5) $ rat))

plot(k3, re3, type="l", xlab = quote(k),  ylim = range(0, re3),
     main = substitute("rel.Err of " ~ (N/Z)^k, list(N=rat3[1], Z=rat3[2])))
if(any((px <- p3^k3) < 2^-1022)) abline(v = min(k3[px < 2^-1022]), lty=3)
abline(h=(-2:2)*2^-53, lty=3-ll, lwd=ll, col=adjustcolor(ll, 1/2))

##===> all these are perfect on Linux Fedora 36 & 38;
## but increasingly (p1, p2, p3)  worse on Windows

## even more extreme
e2 <- -(10:21)
 ps <- 1 - 2 ^ e2
(pM <- 1 - mpfr(2, 256) ^ e2) # 12 'mpfr' numbers of precision  256   bits
## [1]           0.9990234375          0.99951171875          0.999755859375         0.9998779296875
## [5]       0.99993896484375      0.999969482421875      0.9999847412109375     0.99999237060546875
## [9]   0.999996185302734375  0.9999980926513671875  0.99999904632568359375 0.999999523162841796875
stopifnot( (1 - ps) -1 == -ps, ps == pM) ## <==> all ps are exactly representable 1 - 2^{-m}
print(pM, scientific = FALSE, drop0trailing = TRUE)


## The range of interesting k -- such that
## ps^k  does *not* underflow (in double prec) but *is* ~= 10^{-300}:
(ks <- as.integer(round( log(1e-300)/ log(ps) )))
## [1]    707009    1414363    2829071    5658488   11317321  22634987  45270320
## [8]  90540985  181082315  362164975  724330295 1448660934
(k0 <- as.integer(signif(ks, 4)))

cbind(psN <- setNames(seq_along(ps), paste0("1-2^",e2)))
reL  <- lapply(psN, function(i) relEP(ps[i], k0[i] + 0:999, precBits = 2^8))
if(doExtras) { ## much higher precision -- just to show it's *not* needed
    print(system.time(
        reLL <- lapply(psN, function(i) relEP(ps[i], k0[i] + 0:999, precBits = 2^ 14))
    )) # 2.407 sec
    stopifnot(all.equal(reL, reLL, tolerance = 0)) # exactly the same
}

absreL <- lapply(reL, abs)
t(sapply(reL,    summary))
t(sapply(absreL, summary)) # very nice on Linux:
##                 Min.      1st Qu.       Median         Mean      3rd Qu.         Max.
## 1-2^-10 4.698362e-20 1.991498e-17 4.117806e-17 4.292841e-17 6.362789e-17 1.085069e-16
## 1-2^-11 2.644379e-20 2.143401e-17 4.148793e-17 4.352628e-17 6.504935e-17 1.048650e-16
## 1-2^-12 5.185148e-21 2.465148e-17 4.910034e-17 4.719418e-17 6.823981e-17 1.001215e-16
## 1-2^-13 1.732619e-19 2.110998e-17 4.188601e-17 4.175295e-17 6.264898e-17 8.759273e-17
## 1-2^-14 1.028253e-19 2.692311e-17 5.463435e-17 5.234885e-17 7.736664e-17 1.027465e-16
## 1-2^-15 2.173905e-20 1.919186e-17 3.737047e-17 3.683430e-17 5.442273e-17 7.258842e-17
## 1-2^-16 7.547357e-20 1.964997e-17 4.034049e-17 4.078038e-17 6.051778e-17 8.354223e-17
## 1-2^-17 6.468787e-20 2.113894e-17 4.282366e-17 4.184208e-17 6.196293e-17 8.264685e-17
## 1-2^-18 3.053764e-20 2.343907e-17 4.581431e-17 4.515700e-17 6.764085e-17 8.893000e-17
## 1-2^-19 3.693339e-20 1.922513e-17 4.318434e-17 4.283787e-17 6.588232e-17 8.900600e-17
## 1-2^-20 5.398257e-20 2.190955e-17 4.092451e-17 4.095826e-17 6.083059e-17 8.062014e-17
## 1-2^-21 4.708387e-20 2.429069e-17 4.717160e-17 4.758539e-17 7.079720e-17 9.742768e-17

m.absrel <- data.matrix(as.data.frame(absreL, optional=TRUE))
matplot(0:999, m.absrel, type="l", log="y", ylim = c(1e-17, max(m.absrel)))
abline(h = c(1,2,4)*2^-53, lty=2, lwd=2, col=adjustcolor(2, 1/2))

if(.Platform$OS.type == "unix" && !noLdbl) {
    stopifnot(m.absrel < 2.23e-16) # see max(.) == 1.085e-16 < 2^-53 = 1.11022e-16
}

