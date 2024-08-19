## Get  dnt  from 'DPQ' and auto-produce
## i.e., "mpfr-ize" = mpfrize, see also ./dpq-mpfrized.R

if(FALSE) { ## unneeded; dntJKBf() works "alone" already
## dntJKBm() := a version for Rmpfr -- by replacing the final "return line":
nB <- length(bdy <- body(dntJKBm <- dntJKBf))
new <- FALSE # to prevent FP in codetools checks
body(dntJKBm)[[nB]] <- substitute(new("mpfr", THIS), list(THIS = bdy[[nB]]))
rm(nB, bdy, new)
}
## hence deprecated:
dntJKBm <- function (x, df, ncp, log = FALSE, M = 1000) {
    .Deprecated("DPQ::dntJKBf")
    dntJKBf(x, df, ncp, log = FALSE, M = 1000)
}

## Mostly copy_n_paste from ../../DPQ/R/t-nonc-fn.R
##                          ~~~~~~~~~~~~~~~~~~~~~~~
## Orig: ~/R/MM/NUMERICS/dpq-functions/noncentral-t-density-approx_WV.R
## From: Wolfgang Viechtbauer <wviechtb@s.psych.uiuc.edu>
##
##' mpfr-ized version
dtWVm <- function(x, df, ncp=0, log=FALSE) {
    ## even if all (x, df, ncp) are numeric, we use 'mpfr'
    stopifnot(requireNamespace("Rmpfr"))
    if(!any_mpfr(x, df, ncp)) x <- as(x, "mpfr")
    ## we import  mpfr <- Rmpfr::mpfr ; getPrec <- Rmpfr::getPrec ; Const <- Rmpfr::Const
    prec <- max(53L, getPrec(x), getPrec(df), getPrec(ncp))
    pi <- Const("pi", prec = max(64, prec))
    if(!inherits(x,  "mpfr")) x   <- mpfr(x,  prec)
    if(!inherits(df, "mpfr")) df  <- mpfr(df, prec)
    if(!inherits(ncp,"mpfr")) ncp <- mpfr(ncp,prec)
    ln2 <- log(mpfr(2, prec))
    dfx2 <- df + x^2 # = 'f+t^2' in Resnikoff+L.(1957), p.1 (by MM)
    y <- -ncp*x/sqrt(dfx2) # = 'y' in R.+L., p.1
    ## MM(FIXME): cancellation for y >> df  here :
    a <- (-y + sqrt(y^2 + 4*df)) / 2 # NB a = 't' in R.+L., p.25
    dfa2 <- df+a^2 ## << MM(2)
    if(log) {
        lHhmy <- df*log(a) + -0.5*(a+y)^2 +
            0.5*log(2*pi*a^2/dfa2) +
            log1p( - 3*df/(4*dfa2^2) + 5*df^2/(6*dfa2^3))
        lHhmy - (((df-1)/2)*ln2 + lgamma(df/2) + .5*log(pi*df)) +
            -0.5*df*ncp^2/dfx2 + ((df+1)/2)*log(df/dfx2)
    } else { ## MM: cancelled 1/f! = 1/gamma(df+1) in Hh_f(y) =: Hhmy : formula p.25
        Hhmy <- a^df * exp(-0.5*(a+y)^2) *
            sqrt(2*pi*a^2/dfa2) *
            (1 - 3*df/(4*dfa2^2) + 5*df^2/(6*dfa2^3))
        ## formula p.1:  h(f,δ,t) = (....) * Hh_f(-δ t / sqrt(f+t²)) = (....) * Hhmy
        Hhmy / (2^((df-1)/2) * gamma(df/2) * sqrt(pi*df)) *
            exp(-0.5*df*ncp^2/dfx2) * (df/dfx2)^((df+1)/2)
    }
}

###==================== some pbeta()-mpfrized ===================================

## TODO  1) add (and test) R mpfr version R's bpser() in toms708.c ---
##       2) compare with bpser() from package DPQ which should have more flexible variant
##          ~/R/Pkgs/DPQ/R/beta-fns.R  &  ~/R/Pkgs/DPQ/src/bpser.c

### NOTA BENE:  gam1() gamln1() and algdiv() are available as "double precision" ver from DPQ
##  ---------   ~~~~~~ ~~~~~~~~     ~~~~~~~~ {each via C code from toms708.c}
## ==> call the versions here with ...'M' appended

## --->> for now keep "pure R"  versions available  "not installed" (move to documentation ??):
##
##       vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
##  ===> ~/R/Pkgs/DPQ/TODO_R_versions_gam1_etc.R
##       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##
## Here, use Rmpfr-ized versions, calling DPQ's version for double precision argument.

gam1M <- function(a, usePr = NULL) {
    ##    1/Gamma(a+1) - 1   accurately also for |a| << 1

    if(inherits(a, "mpfr")) {
        pr_a <- getPrec(a)
        if(is.null(usePr)) { ## use larger precision
            bits <- ceiling(-log2(abs(a)))
            bits[a == 0] <- 0
            usePr <- pr_a + 2L + pmax(as.integer(bits), 0L)
        } else
            stopifnot(usePr == as.integer(usePr), usePr >= pr_a)
        ## extend precision and "round back" :
        roundMpfr(1/gamma(roundMpfr(a, usePr) + 1) - 1, pr_a)
        ##        ^^^^^^^           ^       ^^^^^^

### __NOT_YET__ FIXME with new DPQ
### } else if(is.numeric(a) && requireNamespace("DPQ")) {
###      DPQ::gam1(a)
    } else {
        warning("not \"mpfr\" -- using direct, possibly inaccurate formula")
                1/gamma(a + 1) - 1
    }
} ## {gam1M}

## gamln1 <- function  ....
lgamma1pM <- function(a, usePr = NULL, DPQmethod = c("lgamma1p", "algam1")) {
##   ---------------------------) {
## == lgamma(a+1) {not clear why this is needed; bpser() used it only for 0 <= a0 < 1
## * ----------------------------------------------------------------------- */
## *     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25 */
## * ----------------------------------------------------------------------- */

### FIXME: in {DPQ}, we have  lgamma1p()
###                           ---------- with Maple derived formulas,
### notably lgamma1p_series() up to 15 terms !! --- should (?) use also with MPFR ?
    if(inherits(a, "mpfr")) {
        ## FIXME -- for really small 'a',
        ##   DPQ::lgamma1p_series() does work with "mpfr"
        pr_a <- getPrec(a)
        if(is.null(usePr)) { ## use larger precision
            bits <- ceiling(-log2(abs(a)))
            bits[a == 0] <- 0
            usePr <- pr_a + 2L + pmax(as.integer(bits), 0L)
        } else
            stopifnot(usePr == as.integer(usePr), usePr >= pr_a)
        ## extend precision and "round back" :
        roundMpfr(lgamma(roundMpfr(a, usePr) + 1), pr_a)
        ##        ^^^^^^^          ^       ^^^^^^
    } else if(is.numeric(a) && requireNamespace("DPQ")) {
        DPQmeth <- get(match.arg(DPQmethod), envir = asNamespace("DPQ"))
        DPQmeth(a)
    } else {
        warning("not numeric, nor \"mpfr\" -- for now using direct, possibly inaccurate formula -- consider DPQ::lgamma1p() , lgamma1p_series()")
        lgamma(a + 1)
    }
}


## NB:   Have  DPQ::algdiv()  double prec.-only which is *based* on TOMS 708 algdiv()
## NB 2: ~/R/Pkgs/DPQ/man/algdiv.Rd -- we mention the relation  logQab(a,b) for we have asymptotic approximations
algdivM <- function(a, b, usePr = NULL) {
## *-------------------------------------------------------------------------*
##      Computation of ln(Gamma(b)/Gamma(a+b))  -- when b >= 8 (in use for bpser() !)
##                                                      ------
## * ----------------------------------------------------------------------- *
    ## MM: ~/R/Pkgs/DPQ/man/algdiv.Rd --> algdiv(a,b) := \log( \Gamma(b) / \Gamma(a+b) )
    ##     ~/R/Pkgs/DPQ/man/lbeta.Rd  --> \log B(a,b) = \log\Gamma(a) - \log Qab ,
    ##     where Qab = Qab(a,b) = \Gamma(a + b) / \Gamma(b)
    ##     ==>    \log Qab(a,b) = \log\Gamma(a + b) - \log\Gamma(b)  =  - algdiv(a,b)
    ##            ^^^^^^^^^^^^^                                         ^^^^^^^^^^^^^
    ##    <==> algdiv(a,b) = lgamma(b)  - lgamma(a+b)   (1)
    ##                     = lbeta(a,b) - lgamma(a)     (2)
    ## Cancellation in (1)  <==> lgamma(a+b) ~= lgamma(b)  <==>  a << b
    ## Cancellation in (2)  <==> lbeta(a,b) ~= lgamma(a)
    ##       <==>  \log[\Gamma(a)\Gamma(b)/\Gamma(a+b)] ~= \log\Gamma(a)
    ##       <==>               \Gamma(b) / \Gamma(a+b) ~= 1
    ##       <==>               \Gamma(a+b) ~= \Gamma(b)
    ##       <==>  a << b

    if((aM <- inherits(a, "mpfr")) |
       (bM <- inherits(b, "mpfr"))) {
        if(!aM) {
            prec <- getPrec(b); a <- mpfr(a, precBits = prec)
        } else if(!bM) {
            prec <- getPrec(a); b <- mpfr(b, precBits = prec)
        } else {
            prec <- max(getPrec(a), getPrec(b))
        }
        if(is.null(usePr)) { ## use "appropriately large" precision
            warning("Choosing #bits for `usePr` is not really tested")
            bits <- ceiling(abs(log2(abs(b)) - log2(abs(a))))
            bits[(a*b) == 0] <- 0
            usePr <- prec + 2L + pmax(as.integer(bits), 0L)
        } else
            stopifnot(usePr == as.integer(usePr), usePr >= prec)
        ## extend precision and "round back" :
        aa <- roundMpfr(a, usePr)
        roundMpfr(lbeta(aa, roundMpfr(b, usePr)) - lgamma(aa), prec)
        ##        ^^^^^^^ ^^          ^        ^^^^^^^^^^^^ ^
    } else if(is.numeric(a) && is.numeric(b) && requireNamespace("DPQ")) {
        DPQ::algdiv(a,b)
    } else {
        warning("not numeric, nor \"mpfr\" -- for now using direct, possibly inaccurate formula -- consider DPQ::logQab_asy()")
        lbeta(a,b) - lgamma(a)
    }
}



## Small subset from ~/R/Pkgs/DPQ/R/dpq-h.R :
##                   ^^^^^^^^^^^^^^^^^^^^^^
.D_0 <- function(log.p) if(log.p) -Inf else 0
.D_1 <- function(log.p) as.integer(!log.p) ## if(log.p) 0 else 1
## Small subset from ~/R/Pkgs/DPQ/R/AAA_consts.R :
##                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^
ML_POSINF <-  Inf
ML_NEGINF <- -Inf
ML_NAN <- NaN

pbeta_ser <- function(q, shape1, shape2, log.p=FALSE,
                      eps = 1e-15, ## R's C: eps = 2. * Rf_d1mach(3); // == DBL_EPSILON ..., but then set to  1e-15
                      errPb = 0, # {0, -1, -2} on input;  {0, 1, 2} on output
                      verbose = FALSE)
{
###___ NB:  ~=  DPQ :: bpser()  which calls the [double prec] C code
###--- ==      --------------   But this one should work with  mpfr-numbers !!

    ## Originally cut'n'paste from R's sources  ~/R/D/r-devel/R/src/nmath/toms708.c, [2022-03-19, lines  512-659]
    ## then  ~/R/Pkgs/DQP/src/bpser.c
    ## -------------------------------------------------------------------------*
    ##   double bpser(double a, double b, double x, double eps, int errPb, int log_p, Rboolean verbose)
    ##  -----------------------------------------------------------------------
    ##  Power SERies expansion for evaluating I_x(a,b) when
    ##         b <= 1 or b*x <= 0.7.   eps is the tolerance used.
    ##  NB: if log.p is TRUE, also use it if   (b < 40  & lambda > 650)   where
    ##                        lambda := a y - b x = (a + b)y - b  = a - (a+b)x   {x + y == 1}
    ##  ----------------------------------------------------------------------- */

    ## such that the remaining code remains closer to the  toms708.c  C code :
    a <- shape1
    b <- shape2
    x <- q
    errPb <- as.integer(errPb) # typically in {-2, -1, 0}

    stopifnot(length(a) == 1, length(b) == 1, length(x) == 1)## TODO:  Vectorize
    if (x == 0.) {
        return(.D_0(log.p))
    }
    # for DPQ: */
    if (x == 1. || a == 0.) return(.D_1(log.p))

# ----------------------------------------------------------------------- */
#             compute the factor  x^a/(a*Beta(a,b)) */
# ----------------------------------------------------------------------- */
    ## double ans, c, a0 = min(a,b);
    a0 <- min(a,b)

    if (a0 >= 1.) { #           ------   1 <= a0 <= b0  ------ */
        z  <- a * log(x) - lbeta(a, b);
        ans  <- if(log.p) z - log(a) else exp(z) / a
    }
    else {
        b0  <- max(a,b)
        if (b0 < 8.) {
            if (b0 <= 1.) { #   ------  a0 < 1  and  a0 <= b0 <= 1  ------ */
                if(log.p) {
                    ans  <- a * log(x);
                } else {
                    ans  <- x^a # pow(x, a);
                    if (ans == 0.) { # once underflow, always underflow .. */
                        if(verbose) cat(sprintf(" bpser(a=%g, b=%g, x=%g): x^a underflows to 0\n",
                                                a,b,x));
                        return(ans)
                    }
                }
                apb  <- a + b;
                if (apb > 1.) {
                    u  <- a + b - 1.;
                    z  <- (gam1M(u) + 1.) / apb;
                } else {
                    z  <- gam1M(apb) + 1.;
                }
                c  <- (gam1M(a) + 1.) * (gam1M(b) + 1.) / z;


                if(log.p) # FIXME ? -- improve quite a bit for c ~= 1 */
                    ans <- ans + log(c * (b / apb))
                else
                    ans <- ans * c * (b / apb)

            } else { #  ------  a0 < 1 < b0 < 8  ------ */

                u  <- lgamma1pM(a0); # was   gamln1(.)

                m  <- as.integer(b0 - 1.) # (int)
                if (m >= 1) {
                    c  <- 1.;
                    for (i in 1:m) {
                        b0 <- b0 - 1
                        c  <- c * b0 / (a0 + b0)
                    }
                    u <- u + log(c);
                }

                z  <- a * log(x) - u;
                b0 <- b0 + -1.; ## => b0 in (0, 7)
                apb  <- a0 + b0;
                if (apb > 1.) {
                    u  <- a0 + b0 - 1.;
                    t  <- (gam1M(u) + 1.) / apb;
                } else {
                    t  <- gam1M(apb) + 1.;
                }

                if(log.p) # FIXME? potential for improving log(t) */
                    ans  <- z + log(a0 / a) + log1p(gam1M(b0)) - log(t)
                else
                    ans  <- exp(z) * (a0 / a) * (gam1M(b0) + 1.) / t
            }

        } else { #              ------  a0 < 1 < 8 <= b0  ------ */

            u  <- lgamma1pM(a0) + algdivM(a0, b0);
            z  <- a * log(x) - u;

            if(log.p)
                ans  <- z + log(a0 / a)
            else
                ans  <- a0 / a * exp(z);
        }
    }
    if(verbose) cat(sprintf(" bpser(a=%g, b=%g, x=%g, log=%d, eps=%g): %s = %.14g;",
                            a,b,x, log.p, eps,
                            if(log.p) "log(x^a/(a*B(a,b)))" else "x^a/(a*B(a,b))", ans))
    if (ans == .D_0(log.p) || (!log.p && a <= eps * 0.1)) {
        if(verbose) cat(" = final answer\n")
        return(ans)
    }
    else if(verbose) cat("\n")

    ## ----------------------------------------------------------------------- */
    ##                      COMPUTE THE SERIES */
    ## ----------------------------------------------------------------------- */
    tol  <- eps / a
    n  <- 0.
    sum  <- 0.
    c  <- 1.
    repeat { ## sum is alternating as long as n < b (<==> 1 - b/n < 0)
        n <- n+1
        c <- c * (0.5 - b / n + 0.5) * x;
        w  <- c / (a + n);
        sum <- sum + w;
        if(!(n < 1e7 && abs(w) > tol)) break
    }
    if(abs(w) > tol) { ## the series did not converge (in time)
        ## warn only when the result seems to matter:
        if(( log.p && !(a*sum > -1. && abs(log1p(a * sum)) < eps*abs(ans))) ||
           (!log.p && abs(a*sum + 1.) != 1.)) {
            if(errPb >= 0) ## caller can specify err_bp = -1  to suppress this warning
                warning(gettextf(" bpser(a=%g, b=%g, x=%g,...) did not converge (n=1e7, |w|/tol=%g > 1; A=%g)",
                        a,b,x, abs(w)/tol, ans), domain=NA)
            errPb  <- 1;
        }
    }
    if(verbose) cat(sprintf("  -> n=%.0f iterations, |w|=%g %s %g=tol:=eps/a ==> a*sum=%g %s -1\n",
                            n, abs(w), if(abs(w) > tol) ">!!>" else "<=", tol,
                            a*sum, if(a*sum > -1.) ">" else "<="))
    if(log.p) {
        if (a*sum > -1.)
            ans <- ans + log1p(a * sum)
        else {
            if(ans > ML_NEGINF) {
                if(errPb >= -1) ## caller can specify err_bp = -2  to suppress both warnings
                    warning(gettextf("pbeta(*, log.p=TRUE) -> bpser(a=%g, b=%g, x=%g,...) underflow to -Inf",
                                     a,b,x), domain=NA)
                errPb  <- 2
            }
            ## FIXME ? rather keep first order term ans = log(x^a/(a*B(a,b))) from above
            ans  <- ML_NEGINF;
        }
    } else if (a*sum > -1.)
        ans <- ans * (a * sum + 1.)
    else ## underflow to
        ans  <- 0.
    ## return
    attr(ans, "err") <- errPb
    ans
} # pbeta_ser {was 'bpser()'} */
## ---------  ------------------------- end cut'n'paste from R sources


## Use Gil_et_al's version of continued fraction ==> ~/R/Pkgs/DPQ/Misc/
##          ~/R/Pkgs/DPQ/R/beta-fns.R  &  ~/R/Pkgs/DPQ/src/bpser.c
pbeta_contFrac <- function(q, shape1, shape2, lower.tail=TRUE, log.p=FALSE) {

    stop("NOT YET!")

}

## "From" ~/R/Pkgs/DPQ/R/t-nonc-fn.R -- simplified, using an MPFR-ized pbeta() from above
##
##' GST23 = Gil et al. (2023) "New asymptotic representations of the noncentral t-distribution"
##' Gil A., Segura J., and Temme N.M. (2023) -- DOI 10.1111/sapm.12609
pntGST23_1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                       j0max = 1e4, # for now
                       pbetaFUN  = pbeta_contFrac,
                       alt = FALSE, # << TODO: find out "automatically" which (T/F) is better
                       verbose = TRUE,
                       ...  ## <-- further arguments passed to pbetaFUN()
                       ) {
    stop("NOT YET!")
}

### ==== Gil et al. (2023) -- DOI 10.1111/sapm.12609 (see above)
###      ----------------
##' Lemma 2, p.861 + Numerics, p.880 (102) --- Integrate Phi(.) === we use our MPFR-Romberg integration !
##'                            ...... ^^^^
if(FALSE) ## __TODO__
pntInorm.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                       integr = c("exp", # trafo e^s = t (102)
                                  "direct"), #  (13)
                       logA = log.p || df >= 286.032165, # ~ empirical boundary
                       ..., # arguments passed to integrate()
                       alt = FALSE,
                       verbose = TRUE)
{
    if(!lower.tail) { ## upper tail .. is *not* optimized for all cases
        if(!log.p)
            1- pntInorm.1(t, df, ncp, lower.tail=TRUE,
                          ..., verbose=verbose)
        else { # log.p: log(1 - F(.))
           if(alt) ## use log scale
               ## alternatively, not equivalently -- see FIXME in log.p case below
               log1mexp(pntInorm.1(t, df, ncp, lower.tail=TRUE, log.p=TRUE,
                                   integr=integr, logA=logA, ..., verbose=verbose))
           else log1p(- pntInorm.1(t, df, ncp, lower.tail=TRUE,
                                   integr=integr, logA=logA, ..., verbose=verbose))
        }
    }

    integr <- match.arg(integr)
    n2 <- df/2 # used "everywhere" below
    ## A_n from (13) .. must use log-scale as soon as n=df is "large"
    ## ===========> FIXME: Use  AnSc() + An.asymp()  from
    ## FIXME               ~/R/Pkgs/DPQ/Misc/pnt-Gil_etal-2023/Giletal23_Ixpq.R
    if(logA)
        lA <- n2*log(n2) - lgamma(n2)
    else
        An <- n2^n2 / gamma(n2)

    switch(integr,
           "exp" = { ## formula (102)
               integrand <- function(....) ...... #_______________________________ FIXME _________________
               I <- integrate(integrand, -Inf, Inf, ...)
           },
           "direct" = { ## formula (13)
               integrand <- function(....) ......
               I <- integrate(integrand, 0, Inf, ...)
           },
           stop("invalid integration method 'integr'=",integr))

    if(verbose) { cat("Integration:\n"); str(I, digits.d = 10) }
    if(log.p) {
        lA + log(I$value)
    } else { ## lower.tail = TRUE,  log.p == FALSE
        if(logA) exp(lA + log(I$value))
        else An  * I$value
    }
}
if(FALSE) ## __TODO__
pntInorm <- Vectorize(pntInorm.1, c("t", "df", "ncp"))


