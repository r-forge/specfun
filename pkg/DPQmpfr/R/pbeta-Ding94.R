## Implement incomplete Beta-distribution according to Ding(1994), Algorithm A, p.452
##
## - take identical arguments as pbeta()
## - non-default (lower.tail, log.p) are *NOT* properly supported
pbetaD94 <- function(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                   eps = 1e-10, itrmax = 100000L, verbose = TRUE) {
    stopifnot(eps >= 0, itrmax <= .Machine$integer.max, # 2^31-1
              ## renaming to (a,b) so notation is more readable:
              length(a <- shape1) == 1, length(b <- shape2) == 1)
    ## FIXME:
    if(!lower.tail) stop("lower.tail=FALSE is not yet supported")
    if(log.p) stop("log.p=TRUE is not supported yet")
    if(inherits(q, "mpfr")) { # ensure compatible precision, and *all* mpfr:
        if(!inherits(a,  "mpfr")) a   <- 0*q + a
        if(!inherits(b,  "mpfr")) b   <- 0*q + b
        if(!inherits(ncp,"mpfr")) ncp <- 0*q + ncp
    }
    ## else they can be double even
    ##
    ## pmin(<M>, *): just in case a guy uses eps=0
    digitsEps <- function(eps, max=1000L) as.integer(pmin(max, round(1 - log10(eps))))
    eps.digits <- digitsEps(eps)
    Fn <- function(x) format(x, digits=eps.digits)

    ## NB: We *could* work in log-scale and actually should  "in extreme cases"
    ## NB(2): FIXME for really large 'ncp' this algorithm "never ends"
    ## Inititialize n=0 :
    ## FIXME: for really large a or b, this overflows to Inf (or Inf/Inf) or underflows to 0
    ##        and we need to work in log-space
    t <- gamma(a+b)/(gamma(a+1)*gamma(b))* q^a * (1-q)^b
    c <- ncp/2 # == \lambda/2
    v <- u <- exp(-c) ## FIXME: if this underflows to 0 we need to work in log-space
    F <- v*t
    n <- 1L
    an <- a+n
    abn <- an+b # = a+b+n
    ab.1 <- abn-1L
    while(an <= abn*q) { ## <- A3
        ## A4:
        u <- u*c/n
        v <- v+u
        t <- t*q*ab.1/an
        F <- F + v*t
        n <- n +1L
        if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in A3-A4 loop")
        an <- an +1L
        ab.1 <- abn
        abn <- abn +1L
    }
    if(verbose) cat("n=",n,"; F=", Fn(F),"  after A3-A4 loop\n", sep="")
    ##
    repeat {
        ## A5 :
        bound <- t*q*ab.1/(an - abn*q)
        ## if(bound <= eps) { # <== as by Ding(1994).. but an absolute bound is *wrong* in the left tail !!
        ## ==> using relative bound :
        if(bound <= eps * F) {
            if(verbose) cat("n=",n," at convergence\n")
            return(F)
        } else if(verbose) cat(sprintf("n=%4d, F=%s, bound= %12s\n", n, Fn(F), format(bound,digits=8)))
        ## A6 (identical to A4, apart from error message):
        u <- u*c/n
        v <- v+u
        t <- t*q*ab.1/an
        F <- F + v*t
        n <- n +1L
        if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in A5-A6 loop")
        an <- an +1L
        ab.1 <- abn
        abn <- abn +1L
    }
    ## never reached
}

### Experimental:
## Quantile of incomplete Beta-distribution according to Ding(1994), Algorithm C, p.453-4
##
## - take identical arguments as qbeta() + extras
## - non-default (lower.tail, log.p) are *NOT* properly supported
dbetaD94 <- function(x, shape1, shape2, ncp = 0, log = FALSE,
                   eps = 1e-10,     # desired accuracy for f()
                   itrmax = 100000L,# maximal number of steps for computing f()
                   verbose = TRUE) {
    stopifnot(eps >= 0, itrmax <= .Machine$integer.max, # 2^31-1
              ## renaming to (a,b) so notation is more readable:
              length(a <- shape1) == 1, length(b <- shape2) == 1)

    if(x == 1) { # (1-x) == 0, not allowed by algorithm
        return(if(b > 1) {
                   if(log) -Inf else 0
               } else if (b == 1) { # f(.) is finite, > 0
                   ## easy : f = sum_{i=0}^Inf dpois(i,lambda/2) * dbeta(1,a+i,1)
                   ##          = a*1 + \sum_{i=0}^Inf i * dpois(i,lambda/2) = a + lambda/2
                   if(log) log(a+ncp/2) else a+ncp/2
               } else { # b < 1
                   Inf
               })
    }
    ## FIXME:
    if(log) stop("log=TRUE is not supported yet")
    if(inherits(x, "mpfr")) { # ensure compatible precision, and *all* mpfr:
        if(!inherits(a,  "mpfr")) a   <- 0*x + a
        if(!inherits(b,  "mpfr")) b   <- 0*x + b
        if(!inherits(ncp,"mpfr")) ncp <- 0*x + ncp
    }
    ## else they can be double even

    ## pmin(<M>, *): just in case a guy uses eps=0
    digitsEps <- function(eps, max=1000L) as.integer(pmin(max, round(1 - log10(eps))))
    Fn <- function(x) format(x, digits=digitsEps(eps))

    ## NB: We *could* work in log-scale and actually should  "in extreme cases"
    ## NB(2): FIXME for really large 'ncp' this algorithm "never ends"
    ## Inititialize n=0 :
    ## FIXME: for really large a or b, this overflows to Inf (or Inf/Inf) or underflows to 0
    ##        and we need to work in log-space
    s <- x^a * (1-x)^b / beta(a,b)
    c <- ncp/2 # == \lambda/2
    v <- u <- exp(-c) ## FIXME: if this underflows to 0 we need to work in log-space
    f <- u*s
    n <- 1L
    ## pre B-3
    an.1 <- a # = a+n-1
    an <- a+1L# = a+n
    abn <- an+b # = a+b+n
    ab.1 <- abn-1L
    while(an <= abn*x) { ## <- B3
        ## B4:
        u <- u*c/n
        v <- v+u
        s <- s*x*ab.1/an.1
        f <- f + u*s
        n <- n +1L
        if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in B3-B4 loop")
        an.1 <- an; an <- an +1L ; ab.1 <- abn ; abn <- abn +1L
    }
    if(verbose) cat("n=",n,"; f=", Fn(f),"  after B3-B4 loop\n", sep="")
    ##
    repeat {
        ## B5 :
        bound <- s*x*ab.1*(1-v)/an.1 # needs 'v' !
        if(bound <= eps * f) {
            if(verbose) cat("n=",n," at convergence\n")
            return(structure(f, iter = n))

        } else if(verbose) cat(sprintf("n=%4d, f=%s, bound= %12s\n",
                                       n, Fn(f), format(bound,digits=8)))
        ## B6 (identical to B4, apart from error message):
        u <- u*c/n
        v <- v+u
        s <- s*x*ab.1/an.1
        f <- f + u*s
        n <- n +1L
        if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in B5-B6 loop")
        if(u*s == 0.) {
            warning("n=",n,": Increments of f underflowed to zero. False convergence")
            break
        }
        an.1 <- an; an <- an +1L ; ab.1 <- abn ; abn <- abn +1L
    }
    ## never reached
} ## {dbetaD94}


## Quantile of incomplete Beta-distribution according to Ding(1994), Algorithm C, p.453-4
##
## - take identical arguments as qbeta() + extras
## - non-default (lower.tail, log.p) are *NOT* properly supported

qbetaD94 <- function(p, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                   delta = 1e-6, # desired accuracy for computing x_p, i.e., the "inverse" via Newton
                   eps = delta^2,  # desired accuracy for computing F() and f()
                   itrmax = 100000L,# maximal number of steps for computing F() and f()
                   iterN = 1000L,    # maximal number of Newton iterations
                   verbose = TRUE) {
    stopifnot(eps >= 0, eps <= delta, itrmax <= .Machine$integer.max, # 2^31-1
              ## renaming to (a,b) so notation is more readable:
              length(a <- shape1) == 1, length(b <- shape2) == 1)
    ## FIXME:
    if(!lower.tail) stop("lower.tail=FALSE is not yet supported")
    if(log.p) stop("log.p=TRUE is not supported yet")
    if(inherits(p, "mpfr")) { # ensure compatible precision, and *all* mpfr:
        if(!inherits(a,  "mpfr")) a   <- 0*p + a
        if(!inherits(b,  "mpfr")) b   <- 0*p + b
        if(!inherits(ncp,"mpfr")) ncp <- 0*p + ncp
    }
    ## else they can be double even

    ## pmin(<M>, *): just in case a guy uses eps=0
    digitsEps <- function(eps, max=1000L) as.integer(pmin(max, round(1 - log10(eps))))
    del.digits <- digitsEps(delta)
    Fn <- function(x) format(x, digits=digitsEps(eps))

    ## NB: We *could* work in log-scale and actually should  "in extreme cases"
    ## NB(2): FIXME for really large 'ncp' this algorithm "never ends"
    ## Inititialize n=0 :
    ## FIXME: for really large a or b, this overflows to Inf (or Inf/Inf) or underflows to 0
    ##        and we need to work in log-space
    cf <- gamma(a+b)/(gamma(a+1)*gamma(b))

    c <- ncp/2 # == \lambda/2
    u0 <- exp(-c) ## FIXME: if this underflows to 0 we need to work in log-space
    x <- 0.5 + 0*p # same type as 'p' (mpfr or double)
    nit <- integer(20L)
    ## Newton Iterations {loop through C4 .. C10} :
    for(jN in seq_len(iterN)) {
        ## C4 :
        n <- 1L
        t <- cf * x^a * (1-x)^b  # (FIXME, overflow / underflow)
        ## NB: We ensure that  1-x > 0  strictly
        s <- a*t/x/(1-x)
        v <- u <- u0
        F <- v*t # CDF
        f <- u*s # PDF
        ## pre-C5:
        an <- a+n
        abn <- an+b # = a+b+n
        ab.1 <- abn-1L
        while(an <= abn*x) { ## C5
            ## C6:
            u <- u*c/n
            v <- v+u
            s <- t*ab.1/(1-x)
            t <- t*x*ab.1/an
            F <- F + v*t
            f <- f + u*s
            n <- n +1L
            if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in C5-C6 loop")
            an <- an +1L
            ab.1 <- abn
            abn <- abn +1L
        }
        if(verbose >= 2)
            cat("n=",n,"; F=", Fn(F)," f=", Fn(f),"  after C5-C6 loop\n", sep="")
        ##
        c1 <- c2 <- FALSE
        repeat {
            ## C7 :
            if(!c1) bound1 <- t* x *  ab.1/(an - abn*x)
            if(!c2) bound2 <- t*(1-v)*ab.1/(1-x)
            ## C8:
            ## Ding(1994) proposed absolute bounds  (bound_j <= eps)  but that's wrong
            ## ==> using relative bounds :
            if((c1 <- bound1 <= eps * F) && (c2 <- bound2 <= eps * f)) {
                if(verbose >= 2) cat("n=",n," at convergence\n")
                break ## --> go to C10
            } else { ## C9 (with parts identical to C6):
                u <- u*c/n
                v <- v+u
                if(c1) {
                    ## !c2 :  F converged; update {s, f} only
                    s <- t*ab.1/(1-x)
                    f <- f + u*s
                    n <- n +1L
                    if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in {C9, f}")
                    ## update an,... and bound2  and go back to C8
                } else if((c2 <- bound2 <= eps * f)) { ## !c1 but c2: f converged; update {t, F} only
                    t <- t*x*ab.1/an
                    F <- F + v*t
                    n <- n +1L
                    if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in {C9, F}")
                    ## update an,... and bound1  and go back to C8
                } else { ## neither F nor f have converged: update all
                    s <- t*ab.1/(1-x)
                    t <- t*x*ab.1/an
                    F <- F + v*t
                    f <- f + u*s
                    n <- n +1L
                    if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in {C9,both} loop")
                }
                if(v*t == 0. && u*s == 0.) {
                    warning("n=",n,": Increments of F and f underflowed to zero. False convergence")
                    break
                }
                an <- an +1L; ab.1 <- abn; abn <- abn +1L
                ## --> update bound1 and/or bound2 and go back to C7 (or C8)
            }
            if(verbose >= 3)
                cat("n=",n,"; F=", Fn(F)," f=", Fn(f),"  inside C7-C9 loop\n", sep="")
        } # {repeat}
        nit[jN] <- n # store #{inner iterations}

        ## C10: Newton step: find new 'x' and check for Newton convergence
        d.x <- (F - p) / f ## {FIXME: overflow / underflow}
        x. <- x - d.x
        if(x. < 0)
            x <- x/2
        else if(1 - x. <= 0) # 1-x *must* not be zero !
            x <- (x+1)/2
        else ## x. := x - d.x is in [0, 1)  {such that 1-x. > 0}
            x <- x.
        if(verbose) {
            usex <- (x. >= 0) && (1 - x. > 0)
            if(verbose >= 2) cat("\n")
            x.string <-
                if(usex) format(x., digits=del.digits)
                else paste0(format(x., digits=5),": not usable; rather x=",
                            format(x, digits = min(8, del.digits)))
            cat(sprintf("Newton d.x=%14.5g; new proposal x=%*s\n",
                        asNumeric(d.x),
                        if(usex) del.digits + 3L else nchar(x.string), x.string))
        }
        if(abs(d.x) <= delta) { ## Newton iterations converged !
            if(verbose) cat("\nNewton converged after", jN, "iterations\n")
            ## and store the "iteration statistics"
            return(structure(x, iterNewton = jN, n.inner = nit[seq_len(jN)]))
        }

    } ## for(jN ..) Newton iterations
}


if(FALSE) {
    ## FIXME: use mpfr

    qbetaD94(0.95, .5, 250, ncp = 0)# already overflows with dbl prec
    ## but using mpfr of course works:
    qbetaD94(1 - 1/mpfr(20 ,64), .5, 250) # 0.007661102468977...
    qbetaD94(1 - 1/mpfr(20,128), .5, 250, eps = 1e-30, delta=1e-20)

    ## Compare with  Table 3  of  Baharev_et_al 2017
    aa <- c(0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 25)
    bb <- c(1:15, 10*c(2:5, 10, 25, 50))

    ## Choose one of these
    p <- 0.95 # = 1 - alpha <<<--- using double precision everywhere below
    p <- 1 - 1/mpfr(20,128) # = 1 - alpha <<<--- using MPFR everywhere below

    delta <- 1e-18
    eps <- delta^2 # default; (delta,eps) are needed for safe-file name
    qbet <- matrix(NA_real_, length(aa), length(bb),
                 dimnames = list(a = formatC(aa), b = formatC(bb)))

    ff <- "/u/maechler/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/Baharev_et_al-2017_table3.txt"
    qbB2017 <- t( data.matrix(read.table(ff)) )
    dimnames(qbB2017) <- dimnames(qbet)

    if(inherits(p, "mpfr")) qbet <- as(qbet, "mpfr")
    print(system.time(
    for(ia in seq_along(aa)) {
        a <- aa[ia]; cat("\na=",a,": b=")
        for(ib in seq_along(bb)) {
            b <- bb[ib]; cat(b," ")
            qbet[ia, ib] <- qbetaD94(p, a, b, ncp = 0, delta=delta)# def.  eps = delta^2
        }
        cat("\n")
    }
    ))# system.time(.)

    ## Safe the expensive computations !!
    if(inherits(p,"mpfr")) {
        prec <- getPrec(p)
        ## FIXME  eps, delta                                    delta eps
        saveF <- paste0("qbetaD94_tab3_Rmpfr_pr", prec,
                        "_",formatC(delta),"_",formatC(eps), ".rds")
        od <- setwd("~/R/MM/Pkg-ex/Rmpfr")
        saveRDS(qbet, file=saveF)
        setwd(od)
    }
    ## number of correct digits:
    summary(nDig <- asNumeric(-log10(abs(1- qbet/qbB2017))))

    matplot(bb, t(asNumeric(-log10(abs(1- qbet/qbB2017)))), type = "o", log="xy")


    ## Looks not good:
    a <-
    ## qbb <- qbetaD94(p, a, b, ncp = 100, delta=1e-18, verbose=3) # verbose=3: too much
    ## but itrmax=100'000 is not sufficient.
    qbb <- qbetaD94(p, a, b, ncp = 100, delta=1e-18, itrmax = 1e8, verbose=2)
    ## Note how silly it is to have a very small 'eps' in a situation were 'x' is still far from the truth
    ## ===> Idea:  Much faster if 'eps' is "large" at the beginning, when the Newton 'd.x' will be inaccurate anyway !!

    pb <- pbetaD94(p, a, b, ncp = 100)


    ## Other ideas:
    ## 1) in case Newton is not usable, be better than  x' = (1+x)/2 {on right hand} : rather use Regula Falsi, or smart unirootR() !
    ## 2) in these cases, use rough estimates, e.g., a few steps of unirootR()    ##    with e.g., eps=1e-3

}

if(FALSE) {
}
