## Implement incomplete Beta-distribution according to Ding(1994), Algorithm A, p.452
##
## - take identical arguments as pbeta()
## - non-default (lower.tail, log.p) are *NOT* properly supported
pbetaD94 <- function(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                     log_scale = (a*b > 0) &&
                         (a+b > 100 || c >= 500), # "arbitrary";  gamma(172) |--> Inf  already
                     ##                                            exp(-750) |--> 0
                     eps = 1e-10, itrmax = 100000L, verbose = FALSE)
{
    stopifnot(eps >= 0, itrmax <= .Machine$integer.max, # 2^31-1
              ## renaming to (a,b) so notation is more readable:
              length(a <- shape1) == 1, length(b <- shape2) == 1)
    ## FIXME:
    if(!lower.tail) stop("lower.tail=FALSE is not yet supported")
    if(log.p) stop("log.p=TRUE is not supported yet")
    useMpfr <- inherits(q, "mpfr")
    formt <- if(useMpfr) Rmpfr::formatMpfr else format
    if(useMpfr) { # ensure compatible precision, and *all* mpfr:
        if(!inherits(a,  "mpfr")) a   <- 0*q + a
        if(!inherits(b,  "mpfr")) b   <- 0*q + b
        if(!inherits(ncp,"mpfr")) ncp <- 0*q + ncp
    }
    c <- ncp/2 # == \lambda/2
    if(verbose) cat("log_scale = ", log_scale, "\n")

    ## else they can be double even
    ##
    ## pmin(<M>, *): just in case a guy uses eps=0
    if(verbose) {
        digitsEps <- function(eps, max=1000L) as.integer(pmin(max, round(1 - log10(eps))))
        eps.digits <- digitsEps(eps)
        Fn <- function(x) formt(x, digits=eps.digits)
    }

    ## NB: We *could* work in log-scale and actually should  "in extreme cases"
    ## NB(2): FIXME for really large 'ncp' this algorithm "never ends"
    ##        and we need to work in log-space <--> 'log_scale'
    ## Inititialize n=0 :
    if(log_scale) { # Work in log scale
        lq <- log(q)
        lt <- lgamma(a+b) - (lgamma(a+1) + lgamma(b)) + a * lq + b *  log1p(-q)
        lv <- lu <- -c
        F <- exp(lv+lt) # F <- v*t
    } else { # (in original scale)
        t <- gamma(a+b)/(gamma(a+1)*gamma(b))* q^a * (1-q)^b
        v <- u <- exp(-c)
        F <- v*t
    }
    n <- 1L
    an <- a+n
    abn <- an+b # = a+b+n
    ab.1 <- abn-1L
    while(an <= abn*q) { ## <- A3
        if(log_scale) { # Work in log scale
            lu <- lu + log(c/n) # u <- u*c/n
            lv <- logspace.add(lv, lu) # lv = log(exp(lv)+exp(lu)) # v <- v+u
            lt <- lt + lq +log(ab.1) - log(an) # t <- t*q*ab.1/an
            F <- F + exp(lv+lt) # F <- F + v*t
        } else { # A4 (in original scale):
            u <- u*c/n
            v <- v+u
            t <- t*q*ab.1/an
            F <- F + v*t
        }
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
        tqa.1 <- if(log_scale) lt + lq + log(ab.1) else t * q * ab.1
        bound <- (if(log_scale) exp(tqa.1) else tqa.1) / (an - abn*q)
        ## if(bound <= eps) { # <== as by Ding(1994).. but an absolute bound is *wrong* in the left tail !!
        ## ==> using relative bound :
        if(bound <= eps * F) {
            if(verbose) cat("n=",n," at convergence\n")
            return(F)
        } else if(verbose >= 2)
            cat(sprintf("n=%4d, F=%s, bound= %12s\n", n, Fn(F), formt(bound, digits=8)))
        ## A6 (identical to A4, apart from error message):
        if(log_scale) { # Work in log scale
            lu <- lu + log(c/n) # u <- u*c/n
            lv <- logspace.add(lv, lu) # lv = log(exp(lv)+exp(lu)) # v <- v+u
            lt <- tqa.1 - log(an) # t <- t*q*ab.1/an
            F <- F + exp(lv+lt) # F <- F + v*t
        } else { # (original scale)
            u <- u*c/n
            v <- v+u
            t <- tqa.1/an
            F <- F + v*t
        }
        n <- n +1L
        if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in A5-A6 loop")
        an <- an +1L
        ab.1 <- abn
        abn <- abn +1L
    }
    ## never reached
} ## {pbetaD94}

## Density of incomplete Beta-distribution according to Ding(1994), Algorithm B, p.453
##
## - take identical arguments as dbeta() + extras
## - non-default (lower.tail, log.p) are *NOT* properly supported
dbetaD94 <- function(x, shape1, shape2, ncp = 0, log = FALSE,
                   eps = 1e-10,     # desired accuracy for f()
                   itrmax = 100000L,# maximal number of steps for computing f()
                   verbose = FALSE)
{
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
    useMpfr <- inherits(x, "mpfr")
    formt <- if(useMpfr) Rmpfr::formatMpfr else format
    if(useMpfr) { # ensure compatible precision, and *all* mpfr:
        if(!inherits(a,  "mpfr")) a   <- 0*x + a
        if(!inherits(b,  "mpfr")) b   <- 0*x + b
        if(!inherits(ncp,"mpfr")) ncp <- 0*x + ncp
    }
    ## else they can be double even

    ## pmin(<M>, *): just in case a guy uses eps=0
    if(verbose) {
        digitsEps <- function(eps, max=1000L) as.integer(pmin(max, round(1 - log10(eps))))
        eps.digits <- digitsEps(eps)
        Fn <- function(x) formt(x, digits=eps.digits)
    }

    ## NB: We *could* work in log-scale and actually should  "in extreme cases"
    ## NB(2): FIXME for really large 'ncp' this algorithm "never ends"
    ## Inititialize n=0 :
    ## FIXME: for really large a or b, this overflows to Inf (or Inf/Inf) or underflows to 0
    ##        and we need to work in log-space  <==> see 'log_scale' in pbetaD94
    s <- x^a * (1-x)^b / beta(a,b)
    c <- ncp/2 # == \lambda/2
    v <- u <- exp(-c) ## FIXME: if this underflows to 0 we need to work in log-space
    f <- u*s  # = 'PDF' in paper
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

        } else if(verbose >= 2) cat(sprintf("n=%4d, f=%s, bound= %12s\n",
                                            n, Fn(f), formt(bound,digits=8)))
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
                   log_scale = (a*b > 0) &&
                       (a+b > 100 || c >= 500), # "arbitrary";  gamma(172) |--> Inf  already
                   ##                                            exp(-750) |--> 0
                   delta = 1e-6, # desired accuracy for computing x_p, i.e., the "inverse" via Newton
                   eps = delta^2,  # desired accuracy for computing F() and f()
                   itrmax = 100000L,# maximal number of steps for computing F() and f()
                   iterN = 1000L,    # maximal number of Newton iterations
                   verbose = FALSE)
{
    stopifnot(eps >= 0, eps <= delta, itrmax <= .Machine$integer.max, # 2^31-1
              ## renaming to (a,b) so notation is more readable:
              length(a <- shape1) == 1, length(b <- shape2) == 1, length(ncp) == 1,
              a >= 0, b >= 0, ncp >= 0)
    ## FIXME:
    if(!lower.tail) stop("lower.tail=FALSE is not yet supported")
    if(log.p) stop("log.p=TRUE is not supported yet")
    useMpfr <- inherits(p, "mpfr")
    formt <- if(useMpfr) Rmpfr::formatMpfr else format
    if(useMpfr) { # ensure compatible precision, and *all* mpfr:
        requireNamespace("Rmpfr", quietly=TRUE)
        if(!inherits(a,  "mpfr")) a   <- 0*p + a
        if(!inherits(b,  "mpfr")) b   <- 0*p + b
        if(!inherits(ncp,"mpfr")) ncp <- 0*p + ncp
    }
    ## else they can be double even

    ## pmin(<M>, *): just in case a guy uses eps=0
    if(verbose) {
        digitsEps <- function(eps, max=1000L) as.integer(pmin(max, round(1 - log10(eps))))
        del.digits <- digitsEps(delta)
        eps.digits <- digitsEps(eps)
        cat(sprintf("del.digits  = digitsEps(eps=delta=%g) =%3d\n", delta, del.digits), sep="",
            sprintf("Fn() digits = digitsEps(eps = %g)     =%3d\n", eps,   eps.digits))
        Fn <- function(x) formt(x, digits=eps.digits)
    }

    ## NB: We *could* work in log-scale and actually should  "in extreme cases"
    ## NB(2): FIXME for really large 'ncp' this algorithm "never ends"

    ## Inititialize n=0 :
    ## For really large a or b, gamma() overflows to Inf (or Inf/Inf) or underflows to 0
    ##        and we need to work in log-space
    c <- ncp/2 # == \lambda/2
    if(verbose) cat("log_scale = ", log_scale, "\n")
    if(log_scale) { # Work in log scale
        lcf <- lgamma(a+b) - (lgamma(a+1) + lgamma(b))
        lu0 <- -c
    } else { # (in original scale)
        cf <- gamma(a+b)/(gamma(a+1)*gamma(b))
        u0 <- exp(-c)
    }
    x <- 0.5 + 0*p # same type as 'p' (mpfr or double)
    nit <- integer(32L) # 32: small number, typically sufficient
    ## Newton Iterations {loop through C4 .. C10} :
    for(jN in seq_len(iterN)) {
        ## C4 :
        n <- 1L
        if(log_scale) { # Work in log scale
            lx <- log(x)
            l1_x <- log1p(-x)
            lt <- lcf + a*lx + b*l1_x ## lt = log(t <- cf * x^a * (1-x)^b)
            ## NB: We ensure that  1-x > 0  strictly
            ls <- log(a) + lt - lx - l1_x ## ls := log(s <- a*t/x/(1-x))
            lv <- lu <- lu0 ## v <- u <- u0
            F <- exp(lv + lt)## F <- v*t # CDF
            f <- exp(lu + ls)## f <- u*s # PDF
        } else { # (in original scale)
            t <- cf * x^a * (1-x)^b
            ## NB: We ensure that  1-x > 0  strictly
            s <- a*t/x/(1-x)
            v <- u <- u0
            F <- v*t # CDF
            f <- u*s # PDF
        }
        ## pre-C5:
        an <- a+n
        abn <- an+b # = a+b+n
        ab.1 <- abn-1L
        while(an <= abn*x) { ## C5
            ## C6:
            if(log_scale) { # Work in log scale
                lu <- lu + log(c/n) # u <- u*c/n
                lv <- logspace.add(lv, lu) # lv = log(exp(lv)+exp(lu)) # v <- v+u
                ls <- lt    +  log(ab.1) - l1_x    # s <- t*ab.1/(1-x)
                lt <- lt + lx +log(ab.1) - log(an) # t <- t*x*ab.1/an
                F <- F + exp(lv+lt) # F <- F + v*t
                f <- f + exp(lu+ls) # f <- f + u*s
            } else { # (in original scale)
                u <- u*c/n
                v <- v+u
                s <- t*ab.1/(1-x)
                t <- t*x*ab.1/an
                F <- F + v*t
                f <- f + u*s
            }
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
            if(log_scale) t <- exp(lt) ## FIXME?  Compute bound *all* in log_space ?
            if(!c1)
                bound1 <- t* x *  ab.1/(an - abn*x)
            if(!c2)
                bound2 <- if(log_scale) t*(-expm1(lv))*ab.1/(1-x)
                          else          t*  (1-v)     *ab.1/(1-x)
            ## C8:
            ## Ding(1994) proposed absolute bounds  (bound_j <= eps)  but that's wrong
            ## ==> using relative bounds :
            if((c1 <- bound1 <= eps * F) && (c2 <- bound2 <= eps * f)) {
                if(verbose >= 2) cat("n=",n," at convergence\n")
                break ## --> go to C10
            } else if(log_scale) { ## C9 (with parts identical to C6) -- Work in log scale
                lu <- lu + log(c/n) # u <- u*c/n
                lv <- logspace.add(lv, lu) # lv = log(exp(lv)+exp(lu)) # v <- v+u
                if(c1) {
                    ## !c2 :  F converged; update {s, f} only
                    ls <- lt + log(ab.1) - l1_x # s <- t*ab.1/(1-x)
                    f <- f + exp(lu+ls)         # f <- f + u*s
                    n <- n +1L
                    if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in {C9, f}")
                    ## update an,... and bound2  and go back to C8
                } else if((c2 <- bound2 <= eps * f)) { ## !c1 but c2: f converged; update {t, F} only
                    lt <- lt + lx +log(ab.1) - log(an) # t <- t*x*ab.1/an
                    F <- F + exp(lv+lt) # F <- F + v*t
                    n <- n +1L
                    if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in {C9, F}")
                    ## update an,... and bound1  and go back to C8
                } else { ## neither F nor f have converged: update all
                    ls <- lt + log(ab.1) - l1_x        # s <- t*ab.1/(1-x)
                    lt <- lt + lx +log(ab.1) - log(an) # t <- t*x*ab.1/an
                    F <- F + exp(lv+lt) # F <- F + v*t
                    f <- f + exp(lu+ls) # f <- f + u*s
                    n <- n +1L
                    if(n > itrmax) stop("n > itrmax=",itrmax," iterations reached in {C9,both} loop")
                }
            } else { ## C9 (with parts identical to C6) -- in original scale
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
            } ## (end C9, original scale)
            an <- an +1L; ab.1 <- abn; abn <- abn +1L
            ## --> update bound1 and/or bound2 and go back to C7 (or C8)
            if(verbose >= 3) {
                cat("n=",n,"; F=", Fn(F)," f=", Fn(f),"  inside C7-C9 loop\n", sep="")
                ## if(verbose >= 4)
                ## convergence if((c1 <- bound1 <= eps * F) && (c2 <- bound2 <= eps * f))
                cat(sprintf(" conv? : (bound1 = %g) %s %g; (bound2 = %g) %s %g\n",
                            bound1, if(c1) "<=" else ">", eps*F,
                            bound2, if(c2) "<=" else ">", eps*f))
            }
        } # {repeat}

        nit[jN] <- n # store #{inner iterations}
        ## C10: Newton step: find new 'x' and check for Newton convergence
        d.x <- (F - p) / f ## {FIXME: overflow / underflow}
        x. <- x - d.x
        x <- if(x. < 0)
                 x/2
             else if(1 - x. <= 0) # 1-x *must* not be zero !
                 (x+1)/2
             else ## x. := x - d.x is in [0, 1)  {such that 1-x. > 0}
                 x.
        if(verbose >= 2) {
            usex <- (x. >= 0) && (1 - x. > 0)
            ## if(verbose >= 2) cat("\n")
            x.string <- ## FIXME: this is not ok for "mpfr"; it becomes
                ##         -----  "<S4 class ‘mpfr1’ [package “Rmpfr”] with 4 slots>"
                if(usex) formt(x., digits=del.digits)
                else paste0(formt(x., digits=5),": not usable; rather x=",
                            formt(x , digits= min(8, del.digits)))
            cat(sprintf("Newton d.x=%14.5g; new proposal x=%*s\n",
                        asNumeric(d.x),
                        if(usex) del.digits + 3L else nchar(x.string), x.string))
        }
        if(abs(d.x) <= delta) { ## Newton iterations converged !
            if(verbose) cat("Newton converged after", jN, "iterations\n")
            ## and store the "iteration statistics"
            return(structure(x, conv = TRUE, iterNewton = jN, n.inner = nit[seq_len(jN)]))
        }
    } ## for(jN ..) Newton iterations
    ## if we "land" here, we did *NOT CONVERGE* ..
    warning(gettextf("qbetaD94() did not converge in %d Newton iterations !", iterN),
            domain=NA)
    structure(x, conv = FALSE, iterNewton = jN, n.inner = nit[seq_len(jN)])
}

## Examples, checks, etc: now in  ../man/pbetaD94.Rd
## and notably in		  ../tests/beta-Ding94.R
