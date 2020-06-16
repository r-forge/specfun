require(Rmpfr)
require(DPQmpfr)

stopifnot(all.equal(DPQ:::dntJKBf(mpfr(0, 64),  5,10), ## gave NaN
                    3.66083172640611114864e-23, tol=1e-20))
## whereas R's 2018-08
dt(0, 5, 10) # is a factor of 2 too large: 7.321663e-23

(dt5.10m <- dntJKBm(mpfr(-4:4, 256), 5, 10))##--> heureka!

stopifnot(all.equal(dt5.10m,  tolerance = 1e-7 ,
                    c(2.604112e-29, 1.40239e-28, 1.423497e-27, 5.449437e-26,
                      3.660832e-23,
                      1.035962e-16, 2.85854e-10, 3.656286e-06, 0.0006233252)))

cbind(x = -4:4,
      dt.log   = dt(-4:4, 5, 10, log=TRUE), # x <= 0: not good!
      log.dt.M = asNumeric(log(dt5.10m)))

## FIXME: ==> more in ~/R/MM/NUMERICS/dpq-functions/dt-ex.R


