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
    mpfr <- Rmpfr::mpfr ; getPrec <- Rmpfr::getPrec
    prec <- max(53L, getPrec(x), getPrec(df), getPrec(ncp))
    pi <- Rmpfr::Const("pi", prec = max(64, prec))
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
