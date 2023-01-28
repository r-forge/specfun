### See also  Rmpfr package containing *many* : ~/R/Pkgs/Rmpfr/R/special-fun.R
##    matches for "[_.A-Za-z0-9]+ <- function" in buffer: special-fun.R
##
##  26: pnorm  <- function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
##  66: dnorm  <- function (x, mean = 0, sd = 1, log = FALSE) {
##  83: dt     <- function (x, df, ncp, log = FALSE) {
## 115: dpois  <- function (x, lambda, log = FALSE,
## 146: dbinom <- function(x, size, prob, log = FALSE,
## 177: dnbinom<- function (x, size, prob, mu, log = FALSE, useLog = any(x > 1e6)) {
## 222: dgamma <- function(x, shape, rate = 1, scale = 1/rate, log = FALSE)
## 457: pbetaI <- function(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,

##' ldexp(f,E) := f * 2^E  (fast and "exact") -- using  DPQ:ldexp  or  Rmpfr::ldexpMpfr
ldexp <- function(f, E) {
    if(is(f, "mpfr"))
        if(packageVersion("Rmpfr") <= "0.9.0") ## workaround
            mpfr(ldexpMpfr(f, E)) else ldexpMpfr(f, E)
    else DPQ::ldexp(f, E)
}

### pnormL* and pnormU*()  from DPQ >= 0.4-2  (2020-10)

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

## These need no code adaption (but use *our* dnorm(), log1mexp, .. ==> environment(.))
pnormU_S53 <- DPQ::pnormU_S53	; environment(pnormU_S53) <- environment()
pnormAsymp <- DPQ::pnormAsymp	; environment(pnormAsymp) <- environment()

### in the future
qnormAsymp <- DPQ::qnormAsymp	; environment(qnormAsymp) <- environment()
if(FALSE) {
## Need to replace  log(M_2PI * x2):
lB <- as.list(body(qnormAsymp))
lB <- as.list(body(qnormAsymp))
hasIf <- vapply(lB, function(o) is.call(o) && o[[1]] == quote(`if`), NA)
## logi [1:10] FALSE FALSE TRUE FALSE FALSE TRUE FALSE ..
## the 2nd if(..) is the one where all the  log(M_2PI * x2) happen
iIf <- which(hasIf)[[2L]]
if(interactive()) {
    str(lB)
    print(hasIf)
    print( lB[[iIf]] )
## if (ord >= 1L) {
"FIXME:"
    ## Here, I need to *add* a line , then all should work fine
    ## >>>- better: Make this an optional argument to qnormAsymp() --> change in DPQ !
'
   M_2PI <- 2 * Const("pi", prec=getPrec(x2))
'
##     x2 <- s2 - log(M_2PI * x2)
##     if (ord >= 2L) {
##         x2 = s2 - log(M_2PI * x2) - 2/(2 + x2)
##         if (ord >= 3L) {
##             x2 = s2 - log(M_2PI * x2) + 2 * log1p(-1/(2 + x2) *
##                 (1 - 1/(4 + x2)))
##             if (ord >= 4L) {
##                 x2 = s2 - log(M_2PI * x2) + 2 * log1p(-1/(2 +
##                   x2) * (1 - 1/(4 + x2) * (1 - 5/(6 + x2))))
##                 if (ord >= 5L) {
##                   x2 = s2 - log(M_2PI * x2) + 2 * log1p(-1/(2 +
##                     x2) * (1 - 1/(4 + x2) * (1 - 1/(6 + x2) *
##                     (5 - 9/(8 + x2)))))
##                 }
##             }
##         }
##     }
## }

    print( lB[[iIf]][[3]] )
    print( lB[[iIf]][[c(3,2)]] )     # log(..)
    print( lB[[iIf]][[c(3,3,3)]] )
}
}#------------------------------------------------
## For now, a copy from (end of) ~/R/Pkgs/DPQ/R/norm_f.R
qnormAsymp <- function(p, lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
                           # ~= log(1-p) -- independent of lower.tail, log.p
                       order,
                       ## begin{added}------------
                       M_2PI = 2* (if(!is.numeric(lp) && is(lp, "mpfr"))
                                       Rmpfr::Const("pi", prec = max(Rmpfr::.getPrec(lp))) else pi),
                       ## end{added}---------------
                       lower.tail = TRUE, log.p = missing(p))
{
    stopifnot(length(order) == 1L, order == (ord <- as.integer(order)), 0L <= ord, ord <= 5L)
    if(missing(p)) { ## log.p unused;  lp = log(1-p)  <==>  e^lp = 1-p  <==>  p = 1 - e^lp
        p. <- -expm1(lp)
        ## swap p <--> 1-p -- so we are where approximation is better
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector
    } else {
        p. <- .D_qIv(p, log.p)
        ## swap p <--> 1-p   (as above)
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector
        p[swap] <- if(log.p) log1mexp(-p[swap]) else 1 - p[swap]
    }
    iFin <- is.finite(lp)
    ##  r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 )
    x2 <- s2 <- -ldexp(lp[iFin], 1) ## = -2*lp = 2s =: xs_0
    if(ord >= 1L) {
        x2 <- s2 - log(M_2PI * x2); ## = xs_1
        if(ord >= 2L) { ## need for (r < 36000.)  <==> s < 36000^2
            x2 = s2 - log(M_2PI * x2) - 2./(2. + x2); ## == xs_2
            if(ord >= 3L) { ## need for (r < 840)
                x2 = s2 - log(M_2PI * x2) + 2*log1p(- 1./(2. + x2)*(1 - 1/(4 + x2))); ## == xs_3
                if(ord >= 4L) { ## need for (r < 109)
                    x2 = s2 - log(M_2PI * x2) + 2*log1p(- 1./(2. + x2)*(1 - 1/(4. + x2)*(1 - 5/(6 + x2)))); ## == xs_4
                    if(ord >= 5L) { ## need for (r < 55)
                        x2 = s2 - log(M_2PI * x2) + 2*log1p( - 1./(2. + x2)*
                                                            (1 - 1/(4. + x2)*
                                                             (1 - 1/(6. + x2)*
                                                              (5 - 9/(8. + x2))))); ## == xs_5
                    }
                }
            }
        }
    }
    R <- -lp # incl. Inf
    R[iFin] <- sqrt(x2)
    R[swap] <- -R[swap]
    ## R[p. == 1/2] <- 0
    R
}


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




