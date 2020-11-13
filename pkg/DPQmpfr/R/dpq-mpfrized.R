## pnormL* and pnormU*()  from DPQ >= 0.4-2  (2020-10)


## Duembgen's lower bound (10), p.6
## { which is strictly better than Komatu(1955)'s lower bound (3) }
pnormL_LD10 <- function(x, lower.tail=FALSE, log.p=FALSE) {
    stopifnot(x > 0)
    ## non-log, upper tail :
    ## 1-Phi(x) >  ~=~ pi*dnorm(x) / ((pi-1)*x + sqrt(2*pi + x^2))
    ## log.p=TRUE and upper tail, i.e.  !lower.tail :
    if(!is.numeric(x) && is(x, "mpfr"))
        pi <- Rmpfr::Const("pi", prec = max(Rmpfr::.getPrec(x)))
    r <- dnorm(x, log=TRUE) - log(x) + log(pi / (pi + sqrt(1 + (2*pi/x)/x) -1))
    if(log.p) {
        if(lower.tail) ## log(1 - exp(r)) = log1mexp(-r)
            log1mexp(-r)
        else
            r
    } else {
        if(lower.tail) -expm1(r) else exp(r)
    }
}

## This one needs no code adaption (but must use "our" dnorm(), etc)
pnormU_S53 <- DPQ::pnormU_S53
environment(pnormU_S53) <- environment()


