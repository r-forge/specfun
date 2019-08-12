### dgamma() --  Catherine Loader's code in R ------

## for now: just the plan
## ~/R/D/r-devel/R/src/nmath/dgamma.c :

## R_D__0 : <==>  (if(log) -Inf else 0)
## R_D__1 : <==>  (if(log) 0 else 1) <==>  !log

dgamma.R <- function(x, shape, scale = 1, log)
{

    if (is.na(x) || is.na(shape) || is.na(scale))
        return (x + shape + scale)

    if (shape < 0 || scale <= 0)
        stop("invalid 'shape' or 'scale'")
    if (x < 0) {
	if(log) -Inf else 0
    } else if (shape == 0) { ##/* point mass at 0 */
	if(x == 0) Inf else if(log) -Inf else 0
    } else if (x == 0) {
	if (shape < 1)
            Inf
        else if (shape > 1) {
            if(log) -Inf else 0
        } else
            if(log) -log(scale) else 1 / scale
    } else if (shape < 1) { ## FIXME: not ok for very small x <==> shape/x overflows to Inf
	pr <- dpois_raw(shape, x/scale, log)
	if(log) pr + (if(shape/x == Inf)log(shape)-log(x) else log(shape/x)) else pr*shape/x
    }
    else {##  shape >= 1
        pr  <- dpois_raw(shape-1, x/scale, log)
        if(log) pr - log(scale) else pr/scale
    }
}

## ~/R/D/r-devel/R/src/nmath/dpois.c   --> dpois_raw()
## // called also from dgamma.c, pgamma.c, dnbeta.c, dnbinom.c, dnchisq.c :
dpois_raw <- function(x, lambda, log)
{
    ##  x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
    ##     lambda >= 0
    R_D__0 <- if(log) -Inf else 0
    R_D__1 <- !log #  (if(log) 0 else 1)
    if (lambda == 0) {
        if(x == 0) R_D__1 else R_D__0
    }
    else if(!is.finite(lambda)) {
        R_D__0 ## including for the case where  x = lambda = +Inf
    } else if (x < 0) {
        R_D__0
    } else if (x <= lambda * (DBL_MIN <- .Machine$double.xmin)) {
        .D_exp(-lambda, log)
    } else if (lambda < x * DBL_MIN) {
	if (!is.finite(x)) ## lambda < x = +Inf
	    R_D__0
	else
            .D_exp(-lambda + x*log(lambda) -lgamma(x+1), log)
    } else
        .D_fexp( 2*pi*x, -stirlerr(x) - bd0(x,lambda), log)
    ##                     ~~~~~~~~      ~~~
}


## ~/R/D/r-devel/R/src/nmath/bd0.c     --> bd0()
bd0 <- function(x, np, verbose = getOption("verbose"))
{
    stopifnot(length(x) == 1)
    if(!is.finite(x) || !is.finite(np) || np == 0.0) {
        ## ML_ERR_return_NAN;
        warning("invalid argument values in  (x, np)")
        return(NaN)
    }

    if(abs(x-np) < 0.1*(x+np)) {
        v <- (x-np)/(x+np) # might underflow to 0
        s <- (x-np)*v
	if(abs(s) < .Machine$double.xmin) ## = DBL_MIN
            return(s)
        ej  <- 2*x*v
	v  <- v*v # // "v = v^2"
	for (j in 1:999) { #/* Taylor series; 1000: no infinite loop
                                        # as |v| < .1,  v^2000 is "zero" */
	    ej <- ej* v ##/ = 2 x v^(2j+1)
	    s_ <- s
	    s  <- s+ ej/(2*j+1)
	    if (s == s_) { ##/* last term was effectively 0 */
                if(verbose)
                    cat(sprintf("bd0(%g, %g): T.series w/ %d terms -> bd0=%g\n",
                                x, np, j, s))
		return(s)
	    }
	}
	warning(gettextf(
            "bd0(%g, %g): T.series failed to converge in 1000 it.; s=%g, ej/(2j+1)=%g\n",
            x, np, s, ej/((2*1000)+1)),
            domain=NA)
    }
    ## else #  | x - np |  is not too small ((or the iterations failed !!))
    x*log(x/np) + np-x
} ## {bd0}


## ~/R/D/r-devel/R/src/nmath/stirlerr.c     --> stirlerr()
## C Code :

##  AUTHOR
##    Catherine Loader, catherine@research.bell-labs.com.
##    October 23, 2000.
##
##  Merge in to R:
 ##	Copyright (C) 2000, The R Core Team

##  DESCRIPTION
##
##    Computes the log of the error term in Stirling's formula.
##      For n > 15, uses the series 1/12n - 1/360n^3 + ...
##      For n <=15, integers or half-integers, uses stored values.
##      For other n < 15, uses lgamma directly (don't use this to
##        write lgamma!)
##
## Merge in to R:
## Copyright (C) 2000, The R Core Team
## R has lgammafn, and lgamma is not part of ISO C
##


##  stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
##              = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
##              = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
##
##  see also lgammacor() in ./lgammacor.c  which computes almost the same!



stirlerr <- function(n) {
    stopifnot(length(n) == 1)
    if (n <= 15.0) {
	nn  <- n + n
	if (nn == (n. <- as.integer(nn)))
            sferr_halves[n. + 1L]
        else ## M_LN_SQRT_2PI = ln(sqrt(2*pi)) = 0.918938..
            lgamma(n + 1.) - (n + 0.5)*log(n) + n - log(sqrt(2*pi))
    }
    else {
        S0 <-  0.083333333333333333333       ## 1/12 */
        S1 <- 0.00277777777777777777778     ## 1/360 */
        S2 <- 0.00079365079365079365079365  ## 1/1260 */
        S3 <- 0.000595238095238095238095238 ## 1/1680 */
        S4 <- 0.0008417508417508417508417508## 1/1188 */
        nn  <- n*n
        if (n > 500) (S0-S1/nn)/n
        else if (n > 80) (S0-(S1-S2/nn)/nn)/n
        else if (n > 35) (S0-(S1-(S2-S3/nn)/nn)/nn)/n
        else ## 15 < n <= 35 :
            (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n
    }
}

##  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
## const static double
sferr_halves <- c(
    0.0, ## n=0 - wrong, place holder only */
    0.1534264097200273452913848,  ## 0.5 */
    0.0810614667953272582196702,  ## 1.0 */
    0.0548141210519176538961390,  ## 1.5 */
    0.0413406959554092940938221,  ## 2.0 */
    0.03316287351993628748511048, ## 2.5 */
    0.02767792568499833914878929, ## 3.0 */
    0.02374616365629749597132920, ## 3.5 */
    0.02079067210376509311152277, ## 4.0 */
    0.01848845053267318523077934, ## 4.5 */
    0.01664469118982119216319487, ## 5.0 */
    0.01513497322191737887351255, ## 5.5 */
    0.01387612882307074799874573, ## 6.0 */
    0.01281046524292022692424986, ## 6.5 */
    0.01189670994589177009505572, ## 7.0 */
    0.01110455975820691732662991, ## 7.5 */
    0.010411265261972096497478567, ## 8.0 */
    0.009799416126158803298389475, ## 8.5 */
    0.009255462182712732917728637, ## 9.0 */
    0.008768700134139385462952823, ## 9.5 */
    0.008330563433362871256469318, ## 10.0 */
    0.007934114564314020547248100, ## 10.5 */
    0.007573675487951840794972024, ## 11.0 */
    0.007244554301320383179543912, ## 11.5 */
    0.006942840107209529865664152, ## 12.0 */
    0.006665247032707682442354394, ## 12.5 */
    0.006408994188004207068439631, ## 13.0 */
    0.006171712263039457647532867, ## 13.5 */
    0.005951370112758847735624416, ## 14.0 */
    0.005746216513010115682023589, ## 14.5 */
    0.005554733551962801371038690  ## 15.0 */
)

