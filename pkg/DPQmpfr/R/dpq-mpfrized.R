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

## This one needs no code adaption (but must use "our" dnorm(), .. ==> environment(.))
pnormU_S53 <- DPQ::pnormU_S53
environment(pnormU_S53) <- environment()


## "Stirling approximation of n!" -- error , i.e., computes
## the log of the error term in Stirling's formula, introduced into R's Mathlib
## by Catherine Loader's  improved formula for dbinom(), dnbinom(), dpois() etc

##' [D]irect formula for stirlerr(), notably adapted to be used with high precision mpfr-numbers
stirlerrM <- function(n, minPrec = 128L) {
    if(notNum <- !is.numeric(n)) {
        precB <- if(isM <- inherits(n, "mpfr"))
                     max(minPrec, .getPrec(n))
                 else if(isZQ <- inherits(n, "bigz") || inherits(n, "bigq"))
                     max(minPrec, getPrec(n))
                 else
                     minPrec
        pi <- Const("pi", precB)
    }
    if(notNum && !isM) {
        if(isZQ)
            n <- mpfr(n, precB)
        else ## the object-author "should" provide a method:
            n <- as(n, "mpfr")
    }
    ## direct formula (suffering from cancellation)
    lgamma(n + 1) - (n + 0.5)*log(n) + n - log(2 * pi)/2
}

##' Few term *asymptotic approximation of stirlerr() -- such it works with bigz, bigq, mpfr
stirlerrSer <- function(n, k) {
    stopifnot(1 <= (k <- as.integer(k)), k <= 11)
    useBig <- (!is.numeric(n) &&
               (inherits(n, "mpfr") || inherits(n, "bigz") || inherits(n, "bigq")))
    if(useBig) {
               ## compute "fully accurate" constants
        frac <- as.bigq
        one <- as.bigz(1)
    } else {
        frac <- `/`
        one <- 1
    }
    ## S_k are bigrational .. perfectly work with "mpfr" or biginteger ("bigz") or bigrational ("bigq")
    S0 <- one/12   # 0.08333333....
    S1 <- one/360  # 0.00277777....
    S2 <- one/1260 # 0.0007936507936..
    S3 <- one/1680 # 0.0005952380952..
    S4 <- one/1188 # 0.00084175084175..
    S5 <- frac(691, 360360) # 0.00191752691752691752695
    S6 <- one/156           # 0.00641025641025641025636
    S7 <- frac(3617,122400) # 0.02955065359477124183007
    S8 <- frac(43867,244188)# 0.17964437236883057316493850
    S9 <- frac(174611,125400) #  1.3924322169059011164274315
    S10<- frac(77683, 5796)   # 13.402864044168391994478957
    ## keep in sync with pkg DPQ's  ~/R/Pkgs/DPQ/R/dgamma.R <<<
    if(is.integer(n))
        n <- as.double(n) # such that  n*n  does not overflow
    n2  <- n*n
    switch(k
        , one/(12*n) # 1
        , (S0-S1/n2)/n # 2
        , (S0-(S1- S2/n2)/n2)/n # 3
        , (S0-(S1-(S2- S3/n2)/n2)/n2)/n # 4
        , (S0-(S1-(S2-(S3- S4/n2)/n2)/n2)/n2)/n # 5
        , (S0-(S1-(S2-(S3-(S4- S5/n2)/n2)/n2)/n2)/n2)/n # 6
        , (S0-(S1-(S2-(S3-(S4-(S5- S6/n2)/n2)/n2)/n2)/n2)/n2)/n # 7
        , (S0-(S1-(S2-(S3-(S4-(S5-(S6 -S7/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n # 8
        , (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-S8/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n # 9
        , (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-S9/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n # 10
        , (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-S10/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n # 11
        )
}




