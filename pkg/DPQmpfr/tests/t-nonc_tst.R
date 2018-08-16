require(Rmpfr)
require(DPQmpfr)

stopifnot(all.equal(dnt.1(mpfr(0, 64),  5,10), ## gave NaN
                    3.66083172640611114864e-23, tol=1e-20))

dntM(mpfr(-4:4, 256), 5, 10)##--> heureka!

## FIXME: ==> more in ~/R/MM/NUMERICS/dpq-functions/dt-ex.R


