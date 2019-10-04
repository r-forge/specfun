#### Functions for the Non-central t-distribution
####                   ===============================================
#### Notation etc from Johnson, Kotz and Balakrishnan (1995)
####                   Chapter 31, '5 Distribution Function', p.514 ff
####                   ===============================================
### Ref.:   (( \cite{JohNKB95} [/u/sfs/bib/Assbib.bib] ))
###  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995).
###  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley Series Prob...

#### Used to be part of ./t-nonc-approx.R, see also ./pt-ex.R
####                      ---------------             -------

## for .DT_val():
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")

### TODO:  E.pnt = E [ . ] = delta * sqrt(nu/2)*gamma(.5*(nu-1))/gamma(.5*nu)
###        -----  is very similar to  b_chi (and can use same asymptotic)

b_chi <- function(nu, one.minus = FALSE, c1 = 341, c2 = 1000)
{
    ## Purpose: Compute  E[ chi_nu ] / sqrt(nu)    --- useful for non-central t
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## Author: Martin Mächler, Date:  6 Feb 1999, 11:21;  Jan.2, 2015
    ## REF: Johnson, Kotz,... , p.520, after (31.26a)
    stopifnot(c1 > 0, c2 >= c1)
    b1 <- function(nu) {
        r <- sqrt(2/nu)*gamma(.5*(nu+1))/gamma(.5*nu)
        if(one.minus) 1-r else r
    }
    b2 <- function(nu) {
        logr <- log(2/nu)/2 + lgamma(.5*(nu+1)) - lgamma(.5*nu)
        ## NOTA BENE: should be much better than 'b1()'  iff one.minus and r << 1
        if(one.minus) -expm1(logr) else exp(logr)
    }
    ## Using boolean vars instead of inefficient  ifelse() :
    r <- nu # will be the result
    r[nu==0] <- if(one.minus) 1 else 0
    ## 3 principal regions
    B1 <- 0  < nu & nu <= c1
    B2 <- c1 < nu & nu <= c2
    BL <- c2 < nu  # "L" : Large
    r[B1] <- b1(r[B1])
    r[B2] <- b2(r[B2])
    r[BL] <- b_chiAsymp(nu[BL], one.minus=one.minus)
    ##       ----------
    r
}

b_chiAsymp <- function(nu, order = 2, one.minus = FALSE)
{
  ## Purpose: Compute  E[ chi_nu ] / sqrt(nu)  --- for "LARGE" nu
  ## ----------------------------------------------------------------------
  ## Arguments: nu >=0  (degrees of freedom)
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 11 Mar 1999 (order = 2);  Aug 21 2018 (order in 1..5)

  ## Abramowitz & Stegun, 6.1.47 (p.257) for a= 1/2, b=0 : __ for  z --> oo ___
  ## gamma(z + 1/2) / gamma(z) ~ sqrt(z)*(1 - 1/(8z) + 1/(128 z^2) + O(1/z^3))
  ## i.e. for b_chi), z = nu/2;  b_chi) = sqrt(1/z)* gamma((nu+1)/2) / gamma(nu/2)
  ##  b_chi) = 1 - 1/(8z) + 1/(128 z^2) + O(1/z^3) ~ 1 - 1/8z * (1 - 1/(16z))
  ##
  ## For order >= 3, I've used Maple expansion etc ==> ~/maple/gamma-exp2.txt
  stopifnot(length(order) == 1L, order == as.integer(order), order >= 1)
  r <- 1/(4*nu)
  r <- r * switch(order, # polynomial order
                  1,			# 1
                  1 - r/2,              # 2
                  1 - r/2*(1+ r* 5),	# 3
                  1 - r/2*(1+ r*(5 - r/4* 21)),	# 4
                  1 - r/2*(1+ r*(5 - r/4*(21 + 399*r))),# 5
                  stop("Need 'order <= 5', but order=",order))
  if(one.minus) r else 1-r
}

## Direct log( sqrt(2/nu)*gamma(.5*(nu+1))/gamma(.5*nu) )
lb_chi00 <- function(nu) {
    n2 <- nu/2
    log(gamma(n2 + 0.5)/ gamma(n2) / sqrt(n2))
}
lb_chi0 <- function(nu) {
    n2 <- nu/2
    lgamma(n2 + 0.5) - lgamma(n2) - log(n2)/2
}

##
lb_chiAsymp <- function(nu, order)
{
    ## Purpose: Asymptotic expansion (nu -> Inf) of  log(b_chi(nu))
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: Aug 23 2018.

    stopifnot(length(order) == 1L, order == as.integer(order), order >= 1)
    ## You can derive the first term from
    ## Abramowitz & Stegun, 6.1.47 (p.257) for a= 1/2, b=0 , or probably
    ## using the Stirling formula for log(Gamma(.))
    if(order == 1)
        return(- 1/4/nu )
    r <- 1/(2*nu) # --> 0
    ## I've used Maple expansion etc ==> ~/maple/gamma-exp2.mw or .tex
    ## from *.tex export of Maple, have
    ## cH := r/2*(-1+(2/3+(-16/5+(272/7+(-7936/9+(353792/11+(-22368256/13+1903757312*r^2*(1/15))*r^2)*r^2)*r^2)*r^2)*r^2)*r^2)
    rr <- r*r # = r^2
    O <- rr*0 # in correct (Rmpfr) precision
    -r/2 *
        switch(order, # polynomial order {written to use full prec w Rmpfr}
               1,			      # 1  (degree 1)
               1 - rr*2/3,                   # 2  (deg.   3)
               1 - rr*((O+2)/3 - rr*16/5),   # 3  (deg.   5)
               1 - rr*((O+2)/3 - rr*((O+16)/5 - rr*272/7)),   # 4  (deg. 7)
               1 - rr*((O+2)/3 - rr*((O+16)/5 - rr*(O+272)/7 - rr*7936/9)),# 5 (deg. 9)
               stop("Currently need 'order <= 5', but order=",order))
}


## was called 'pt.appr1'
pntLrg <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE) {
    ## Approximate pt()  [non-central t] --- using the "extreme case" formula
    ##  from C's pnt.c  and  pntR() below:
    ##
    ## if (df > 4e5 || del*del > 2*M_LN2*(-(DBL_MIN_EXP))) {
    ## /*-- 2nd part: if del > sqrt(2*log(2)*1021) = 37.62.., then p=0 below
    ## FIXME: test should depend on `df', `tt' AND `del' ! */
    ## /* Approx. from	 Abramowitz & Stegun 26.7.10 (p.949) */

    ## Vectorizing (in t) {{is this correct? -- FIXME??}}
    n <- max(length(t), length(df), length(ncp))
    if(length( t ) != n) t   <- rep_len(t,  n)
    if(length( df) != n) df  <- rep_len(df, n)
    if(length(ncp) != n) ncp <- rep_len(ncp,n)

    neg <- t < 0
    t  [neg] <- -   t[neg]
    ncp[neg] <- - ncp[neg]
    s <- 1/(4*df)
    pnorm(t*(1 - s), mean = ncp, sd = sqrt(1 + t*t*2*s),
          lower.tail = (lower.tail != neg), log.p=log.p)
}

## was called 'pt.appr'
pntJW39.0 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE)
{
  ## Purpose: Jennett & Welch (1939) approximation to non-central t
  ##          see Johnson, Kotz, Bala... Vol.2, 2nd ed.(1995)
  ##                                     p.520, after (31.26a)
  ## -- still works FAST for huge ncp
  ##    but has *wrong* asymptotic tail (for |t| -> oo)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?pt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date:  6 Feb 1999, 11 Mar 1999

  b <- b_chi(df)
  ##   =====
  ## FIXME:  (1 - b^2) below suffers from severe cancellation!
  pnorm((t*b - ncp)/sqrt(1+ t*t*(1 - b*b)),
        lower.tail = lower.tail, log.p = log.p)
}

pntJW39 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE)
{
  ## Purpose: Jennett & Welch (1939) approximation to non-central t
  ##          see Johnson, Kotz, Bala... Vol.2, 2nd ed.(1995)
  ##                                     p.520, after (31.26a)
  ## -- still works FAST for huge ncp
  ##    but has *wrong* asymptotic tail (for |t| -> oo)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?pt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date: Feb/Mar 1999; 1-b^2 improvement: Aug 2018
  ._1_b <- b_chi(df, one.minus=TRUE)# == 1 - b -- needed for good (1 - b^2)
  ##       =====
  b <- 1 - ._1_b
  ## (1 - b^2) == (1 - b)(1 + b) = ._1_b*(2 - ._1_b)
  pnorm((t*b - ncp)/sqrt(1+ t*t * ._1_b*(1 + b)),
        lower.tail = lower.tail, log.p = log.p)
}

c_dt <- function(nu) {
    ## Purpose: The constant in  dt(t, nu, log=TRUE) =
    ##         = \log f_{\nu}(t) = c_{\nu} - \frac{\nu+1}{2}\log(1 + x^2/\nu)
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  1 Nov 2008, 14:26

    warning("FIXME: current c_dt() is poor -- base it on lb_chi(nu) !")

    ## Limit nu -> Inf: c_{\nu} \to - \log(2\pi)/2 = -0.9189385
    r <- nu
    nu. <- nu[I <- nu < 200]
    r[I] <- log(gamma((nu.+1)/2)/gamma(nu./2)/sqrt(pi*nu.))
    nu. <- nu[I <- !I & nu < 1e4]
    r[I] <- lgamma((nu.+1)/2) - lgamma(nu./2) - log(pi*nu.)/2
    I <- nu >= 1e4
    r[I] <- c_dtAsymp(nu[I])
    r
}

c_dtAsymp <- function(nu)
{
    ## Purpose: Asymptotic of c_dt -- via Abramowitz & Stegun, 6.1.47 (p.257)
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ##
    ## FIXME: This is trivially   -log(2*pi)/2  +  log(b_chi(nu))
    ##        and I've computed good asymptotics for
    ##    lb_chi(nu) := log(b_chi(nu))   above
    ##
    ##
    ## -log(2*pi)/2 -1/(4*nu) * (1 + 5/(96*nu^2))
    ##                               ^^^^^^^^^^^ not quite ok
    warning("this is poor -- use lb_chi(nu) !!")
    -log(2*pi)/2 - 1/(4*nu)
}

c_pt <- function(nu)
{
  ## Purpose: the asymptotic constant in log F_nu(-x) ~= const(nu) - nu * log(x)
  ##          where F_nu(x) == pt(x, nu)
### FIXME == Source / Reference for the above left tail statement?
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 Nov 2008, 09:34
  warning("use better c_dt()") ## ==> use lb_chi(nu) !!
  c_dt(nu) + (nu-1)/2 * log(nu)
}


## MM:  My conclusion of experiments in ../tests/pnt-precision-problem.R :
## -------              --------------------------
## Large t (>0 or < 0)  *MUST* get a new algorithm !
## -------              --------------------------

pntR1  <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                   itrmax = 1000, errmax = 1e-12, verbose = TRUE)
{
    ## Purpose: R version of the series used in pnt() in
    ##          ~/R/D/r-devel/R/src/nmath/pnt.c
    ##
    ## ----------------------------------------------------------------------
    ## Arguments: same as  pt()
    ## Author: Martin Mächler, Date:  3 Apr 1999, 23:43

    stopifnot(length(t) == 1, length(df) == 1, length(ncp) == 1,
              df > 0, ncp >= 0) ## ncp == 0 --> pt()
    ## Make this workable also for "mpfr" objects --> use format() in cat()
    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(t, df, ncp)
    if(!isN) {
	if(verbose) cat("some 'non-numeric arguments .. fine\n")
	if(isMpfr) {
	    stopifnot(requireNamespace("Rmpfr"))
            getPrec <- Rmpfr::getPrec
	    prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
	    pi <- Rmpfr::Const("pi", prec = max(64, prec))
	}
    }

    if (t >= 0) {
        negdel <- FALSE; tt <- t; del <- ncp
    } else {
        if(verbose) cat("t < 0  ==> swap delta := -ncp\n")
	negdel <- TRUE ; tt <- -t; del <- -ncp
    }

    if (df > 4e5 || del*del > 2*log(2)*(-.Machine$double.min.exp)) {
        ## /*-- 2nd part: if del > 37.6403, then p=0 below
        ## FIXME: test should depend on `df', `tt' AND `del' ! */
        ## /* Approx. from Abramowitz & Stegun 26.7.10 (p.949) */
	s <- 1./(4.*df)
        pnt.. <- pnorm(tt*(1 - s), del, sqrt(1. + tt*tt*2.*s),
                       lower.tail = (lower.tail != negdel), log.p=log.p)
        cat("large 'df' or \"large\" 'ncp' -- C code would return pnorm(*) =",
            format(pnt.., digits=16), "\n")
        ### FIXME: see above !!
        return(pnt..)
    }

    ## initialize twin series */
    ## Guenther, J. (1978). Statist. Computn. Simuln. vol.6, 199. */

    x <- t * t
    rxb <- df/ (x + df) ## := (1 - x) {x : below} -- but more accurately
    x   <- x / (x + df) # in [0,1) -- x == 1 for very large t */
    ## (1 - x) = df / (x + df)

    if(verbose)
        cat("pnt(t=",format(t),", df=",format(df),", delta=",format(ncp),") ==> x= ",
            format(x),":")

    tnc <- 0 # tnc will be the result */
    if (x > 0) { ## <==>  t != 0 */
	lambda <- del * del
	p <- 0.5 * exp(-0.5 * lambda)
	if(verbose) cat("\t p=",format(p),"\n")

	if(p == 0.) { # underflow! */

## FIXME: should "analytically  factor-out the small 0.5 * exp(-0.5 * lambda)
## -----  and re-multiply at end --- where we need *logspace* addition of
## in      tnc <-  tnc*0.5 * exp(-0.5 * lambda)  +  pnorm(- del, 0, 1)
## in      tnc <-  exp(log(tnc*0.5 * exp(-0.5 * lambda))  + exp(pnorm(- del, 0, 1, log=TRUE))
## in      log.tnc <-  log(exp(AA) + exp(BB)) === copula:::lsum(rbind(AA, BB))
##            where  AA <- log(0.5) + log(tnc) - lambda/2
##            and    BB <- pnorm(- del, log.p=TRUE)


            ##========== really use an other algorithm for this case !!! */
            cat("\n__underflow__ -- p = 0 -- |ncp| = |delta|  too large\n\n")
            ## |delta| too large */
            return(0)
	}
        if(verbose >= 2)
            cat(  " it   1e5*(godd,  geven)         p          q         s ",
                ## 1.3 1..4..7.9 1..4..7.9|1..4..7.9  1..4..7.9  1..4..7.9_ */
                  "           pnt(*)    D(pnt)     errbd\n", sep=''
                ## 1..4..7..0..3..67 1..4..7.9 1..4..7.9*/
                )

	q <- sqrt(2/pi) * p * del
	a <- 0.5
	b <- 0.5 * df
	s <- 0.5 - p # but that may suffer from cancallation:
        ## s =  0.5 - p = 0.5*(1 - exp(-L)) =  -0.5*expm1(-L))
        if(s < 1e-7)
            s <- -0.5 * expm1(-0.5 * lambda)
	## was: rxb <- (1 - x) ^ b
        ## (1 - x) = df / (x + df)
        rxb <- rxb ^ b
	albeta <- .5*log(pi) + lgamma(b) - lgamma(0.5 + b)
	xodd <- if(isN) pbeta(x, a, b) else pbetaRv1(x, a, b)
        ##                                  -------- ./beta-gamma-etc/pbetaR.R
	godd <- 2 * rxb * exp(a * log(x) - albeta)
	xeven <- if(b*x <= .Machine$double.eps) b*x else 1 - rxb
        ## xeven = 1 - (1 - x)^b = b*x - choose(b,2) * x^2  + O((bx)^3)
	geven <- b * x * rxb
	tnc <- p * xodd + q * xeven
        errbd <- Inf

        ## repeat until convergence or iteration limit */
	for(it in 1:itrmax) {
	    a <- a+1
	    xodd  <- xodd - godd
	    xeven <- xeven - geven
	    godd  <- godd * x * (a + b - 1) / a
	    geven <- geven* x * (a + b - 0.5) / (a + 0.5)
	    p <- p * lambda / (2 * it)
	    q <- q * lambda / (2 * it + 1)
	    tnc <- tnc + p * xodd + q * xeven
	    s <- s - p
	    if(s < -1e-10){## happens e.g. for (t,df,delta)=(40,10,38.5), after 799 it.*/
                cat("IN loop (it=",it,"): precision NOT reached\n")
		## ML_ERROR(ME_PRECISION)
                if(verbose)
                    cat("s =",format(s)," < 0 !!! ---> non-convergence!!\n")
		break # goto finis
	    }
	    if(s <= 0 && it > 1) break # goto finis
	    errbd <- 2 * s * (xodd - godd)
            if(verbose >= 2)
                cat(sprintf(paste("%3d %#9.4g %#9.4g",
                                  "%#10.4g %#10.4g %#9.4g %#17.15g %#9.4g %#9.4g\n",
                                  sep="|"),
                            it, 1e5*godd, 1e5*geven,
                            p,q, s, tnc, p * xodd + q * xeven ,errbd))
	    if(abs(errbd) < errmax) break # convergence
	}
        ## non-convergence:*/
        if(abs(errbd) >= errmax && s > 0)
            cat("end loop: precision NOT reached\n") ## ML_ERROR(ME_PRECISION)
    } ## x > 0
    tnc0 <- tnc
    tnc <- tnc + pnorm(- del, 0, 1)

    ## was -- if (negdel) 1 - tnc else tnc
    lower.tail <- lower.tail != negdel ## xor

    if(verbose)
        cat(sprintf("%stnc{sum} = %.12g, tnc+pnorm(-del) = %.12g, lower.t = %d\n",
                    if(verbose == 1 && x > 0) sprintf("%d iter.: ", it) else "\\--> ",
                    tnc0, tnc, lower.tail))

    if(tnc > 1 - 1e-10 && lower.tail)
        cat("finis: precision NOT reached\n")

    .DT_val(min(tnc, 1.), lower.tail, log.p)
}# pntR1()

pntR  <- Vectorize(pntR1, c("t", "df", "ncp"))
##==

##' Simple vector version of  copula:::lsum() -- ~/R/Pkgs/copula/R/special-func.R
##' Properly compute log(x_1 + .. + x_n) for given log(x_1),..,log(x_n)
##' Here, x_i > 0  for all i
##'
##' @title Properly compute the logarithm of a sum
##' @param lx n-vector of values log(x_1),..,log(x_n)
##' @param l.off the offset to substract and re-add; ideally in the order of
##'        the maximum of each column
##' @return log(x_1 + .. + x_n) = log(sum(x)) = log(sum(exp(log(x))))
##'         = log(exp(log(x_max))*sum(exp(log(x)-log(x_max))))
##'         = log(x_max) + log(sum(exp(log(x)-log(x_max)))))
##'         = lx.max + log(sum(exp(lx-lx.max)))
##' @author Martin Maechler (originally joint with Marius Hofert)
##'
##' NB: If lx == -Inf for all should give -Inf, but gives NaN / if *one* is +Inf, give +Inf
lsum <- function(lx, l.off = max(lx)) {
    if (is.finite(l.off))
        l.off + log(sum(exp(lx - l.off)))
    else if(missing(l.off) || is.na(l.off) || l.off == max(lx))
        l.off
    else stop("'l.off  is infinite but not == max(.)")
}


##' Simple vector version of  copula:::llsum() -- ~/R/Pkgs/copula/R/special-func.R
##' Properly compute log(x_1 + .. + x_n) for given
##' log(|x_1|),.., log(|x_n|) and corresponding signs sign(x_1),.., sign(x_n)
##' Here, x_i is of arbitrary sign
##' @title compute logarithm of a sum with signed large coefficients
##' @param lxabs n-vector of values log(|x_1|),..,log(|x_n|)
##' @param signs corresponding signs sign(x_1), .., sign(x_n)
##' @param l.off the offset to substract and re-add; ideally in the order of max(.)
##' @param strict logical indicating if it should stop on some negative sums
##' @return log(x_1 + .. + x_n)
##'         log(sum(x)) = log(sum(sign(x)*|x|)) = log(sum(sign(x)*exp(log(|x|))))
##'         = log(exp(log(x0))*sum(signs*exp(log(|x|)-log(x0))))
##'         = log(x0) + log(sum(signs* exp(log(|x|)-log(x0))))
##'         = l.off   + log(sum(signs* exp(lxabs -  l.off  )))
##' @author Martin Maechler (originally joint with Marius Hofert)
lssum <- function (lxabs, signs, l.off = max(lxabs), strict = TRUE)
{
    sum. <- sum(signs * exp(lxabs - l.off))
    if (any(is.nan(sum.) || sum. <= 0))
        (if (strict) stop else warning)("lssum found non-positive sums")
    l.off + log(sum.)
}


### Simple inefficient but hopefully correct version of pntP94..()
### This is really a direct implementation of formula
### (31.50), p.532 of  Johnson, Kotz and Balakrishnan (1995)

pnt3150.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE, M = 1000,
                      verbose = TRUE) {
    stopifnot(length(t) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              is.numeric(M), M >= 1, M == round(M))
    if (t < 0)
        return(pnt3150.1(-t, df, ncp = -ncp, lower.tail = !lower.tail,
                         log.p=log.p, M=M))

    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(t, df, ncp)
    if(!isN) {
        if(isMpfr) {
            stopifnot(requireNamespace("Rmpfr"))
            ## if(!exists("pbetaRv1", mode="function"))
            ##     source("~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/pbetaR.R")
            getPrec <- Rmpfr::getPrec
            prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
            pbeta <- Vectorize(pbetaRv1, "shape2") # so pbeta(x, p, <vector q>) works
            if(FALSE) {
                pi <- Rmpfr::Const("pi", prec = max(64, prec))
                dbeta <- function(x, a,b, log=FALSE) {
                    lval <- (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b)
                    if(log) lval else exp(lval)
                }
            }
        }
    }

    ## hence now, t >= 0 :
    x <- df / (df + t^2)
    lambda <- 1/2 * ncp^2
    ## Cheap and stupid -- and not ok for moderately large ncp !
    j <- 0:M
    ## terms <- (ncp/sqrt(2))^j / gamma(j/2 + 1) * pbeta(x, df/2, (j+1)/2)
    ## log(.):
    l.terms <- j*(log(ncp)- log(2)/2) - lgamma(j/2 + 1) + pbeta(x, df/2, (j+1)/2, log.p=TRUE)
    if(any(ina <- is.na(l.terms))) {
        ina <- which(ina)
        if(ina[length(ina)] == M+1 && all(diff(ina) == 1)) {
            ## NaN [underflow] only for j >= j_0
            cat("good: NaN's only for j >= ", ina[1],"\n")
            l.terms <- l.terms[seq_len(ina[1]-1)]
        }
    }
    if(verbose) {
        cat(sprintf(" log(terms)[1:%d] :\n", length(l.terms))); print(summary(l.terms))
        elt <- exp(l.terms)
        cat(sprintf("sum(exp(l.terms)) , sum(\"sort\"(exp(l.terms))) = (%.16g, %.16g)\n",
                    sum(elt), sum(exp(sort(l.terms)))))
        cat(sprintf("log(sum(exp(l.terms))) , lsum(l.terms), rel.Delta = (%.16g, %.16g, %.5e)\n",
                    log(sum(elt)), lsum(l.terms), 1 - log(sum(elt))/lsum(l.terms)))
        cat("exp(-delta^2/2) =", format(exp( - ncp^2/2)),"\n")
    }
    ## P(..) = 1 - exp(-lambda)/2 * Sum == 1 - exp(-LS)
    ## where  LS := lambda + log(2) - log(Sum)
    ##  exp(-lambda)/2 * Sum = exp(-lambda - ln(2) + log(Sum))

    ## For non-small ncp, and small t (=> very small P(.)),
    ##  e.g., pt(3, 5, 10) have huge cancellation here :
    LS <- lambda + log(2) - lsum(l.terms)

    if(log.p) {
        if(lower.tail) ## log(1 - exp(-LS)) = log1mexp(LS)
            log1mexp(LS)
        else ## upper tail: log(1 - (1 - exp(-LS))) = -LS
            -LS
    } else {
        if(lower.tail) -expm1(-LS) else exp(-LS)
    }
}
pnt3150  <- Vectorize(pnt3150.1, c("t", "df", "ncp"))

### New version of pntR1(), pntR() ----  Using the  Posten (1994) algorithm
pntP94.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                     itrmax = 1000, errmax = 1e-12, verbose = TRUE)
{
    stopifnot(length(t) == 1, length(df) == 1, length(ncp) == 1,
              df > 0, ncp >= 0) ## ncp == 0 --> pt()

    Cat <- function(...) if(verbose > 0) cat(...)
    ## Make this workable also for "mpfr" objects --> use format() in cat()
    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(t, df, ncp)
    if(!isN) {
        if(isMpfr) {
            stopifnot(requireNamespace("Rmpfr"))
            ## if(!exists("pbetaRv1", mode="function"))
            ##     source("~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/pbetaR.R")
            getPrec <- Rmpfr::getPrec
            prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
            pi <- Rmpfr::Const("pi", prec = max(64, prec))
            pbeta <- pbetaRv1
            dbeta <- function(x, a,b, log=FALSE) {
                lval <- (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b)
                if(log) lval else exp(lval)
            }
        }
        ## if(FALSE) ## needed for printing mpfr numbers {-> pkg Rmpfr}, e.g.
	## .N <- if(isMpfr) Rmpfr::asNumeric else as.numeric
        Cat("Some 'non-numeric arguments .. fine\n")
    }

    neg.t <- (t < 0)
    if(neg.t) {
        if(verbose) cat("t < 0  ==> swap delta := -ncp\n")
	## tt <- -t
        del <- -ncp
    } else {
        ## tt <- t
        del <- ncp
    }

    x <- t * t
    ## x := df / (df + t^2)   <==>  (1-x) =: rxp = t^2 / (df + t^2)
    x   <- df/ (x + df) ## in (0, 1] :  x == 1 for very large t
    rxb <- x / (x + df) # == (1 - x)  in [0,1)
    ## (1 - x) = df / (t^2 + df)

    Cat("pnt(t=",format(t),", df=",format(df),", delta=",format(ncp),") ==> x= ",
        format(x),":")

    lambda <- del * del / 2
    B  <- pbeta(x, df/2, 1/2)
    BB <- pbeta(x, df/2, 1  )
    x.1x <- x*rxb # ==  x (1 - x)
    S  <- 2*x.1x*dbeta(x, df/2, 1/2)
    SS <-   x.1x*dbeta(x, df/2, 1  )

    L <- lambda ## at the end multiply 'Sum' with exp(-L)
    ## "FIXME" e.g. for large lambda -- correct already: rescale (B,BB) and/or (S,SS)
    D <- 1
    E <- D*del*sqrt(2/pi)

    Sum <- D*B + E * BB # will be the result

    ## repeat until convergence or iteration limit */
    for(i in 1:itrmax) {
        B  <-  B + S
        BB <- BB + SS
        D <- lambda/ i        * D
        E <- lambda/(i + 1/2) * E
        term <- D*B + E * BB
        Sum <- Sum + term
        if(verbose >= 2)
            cat("D=",format(D),", E=",format(E),
                "\nD*B=",format(D*B),", E*BB=",format(E*BB),
                "\n -> term=",format(term),", sum=",format(Sum),
                ";\n rel.chng=", format(term/Sum),"\n")
        if(abs(term) <= errmax * abs(Sum))
            break ## convergence

        if(i < itrmax) {
            i2 <- i*2
            S  <- rxb * (df+i2-1)/(i2+1) * S
            SS <- rxb * (df+i2  )/(i2+2) * SS
        }
        else
            warning("Sum did not converge with ", itrmax, "terms")
    }

    if(neg.t) lower.tail <- !lower.tail

    ## FIXME 1: for large lambda [Sum maybe too large (have overflown) as well!]
    ## FIXME 2: if (log.p) do better anyway!

    ## iP = 1 - P(..) = exp(-L)*Sum/2 == exp(-LS), as
    ## exp(-L) * S/2 = exp(-L)*exp(log(S/2)) = exp(-L + log(S) - log(2)) =: exp(-LS)
    ## where  LS := L +log(2) - log(S)
    LS <- L + log(2) - log(Sum)
    if(verbose)
        cat( if(verbose == 1) sprintf("%d iter.: ", i) else "\\--> ",
            "Sum = ", format(Sum),
            "LS = -log(1-P(..))=", format(LS), "\n")

    if(log.p) {
        if(lower.tail) ## log(1 - exp(-L)* S/2) = log(1 - exp(-LS)) = log1mexp(LS)
            log1mexp(LS)
        else
            -LS
    } else {
        if(lower.tail) -expm1(-LS) else exp(-LS)
    }
}
pntP94 <- Vectorize(pntP94.1, c("t", "df", "ncp"))

### Should implement
pntChShP94.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                   itrmax = 1000, errmax = 1e-12, verbose = TRUE) {
    stop("not yet ...")
    ## ~/save/papers/Numerics/Chattamvelli+Sh_noncentr-t_1994.pdf

    ## Chattamvelli, R. and Shanmugam, R. (1994)
    ## An enhanced algorithm for noncentral t-distribution,
    ## \emph{Journal of Statistical Computation and Simulation} \bold{49}, 77--83.
}

pntChShP94 <- Vectorize(pntChShP94.1, c("t", "df", "ncp"))


###--- dnt() --- for the non-central *density*
###                      =====================
### Wikipedia (or other online sources) have no formula for the density

### Johnson, Kotz and Balakrishnan (1995) [2nd ed.] have
###  (31.15) [p.516] and (31.15'), p.519 -- and they  contradict by a factor  (1 / j!)
###  in the sum ---> see below


if(FALSE) {## Just for historical reference :
           ## =============================
dntRwrong1 <- function(x, df, ncp, log = FALSE, M = 1000, check=FALSE, tol.check = 1e-7)
{
    ## R's source ~/R/D/r-devel/R/src/nmath/dnt.c  claims -- from 2003 till 2014 -- but *WRONGLY*
    ## *	   f(x, df, ncp) =
    ## *		df^(df/2) * exp(-.5*ncp^2) /
    ## *		(sqrt(pi)*gamma(df/2)*(df+x^2)^((df+1)/2)) *
    ## *		sum_{k=0}^Inf  gamma((df + k + df)/2)*ncp^k /
    ## *				prod(1:k)*(2*x^2/(df+x^2))^(k/2)
    stopifnot(length(x) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              ncp >= 0, df > 0,
              is.numeric(M), M >= 1, M == round(M))

    if(check)
     fac <- df^(df/2) * exp(-.5*ncp^2) / (sqrt(pi)*gamma(df/2)*(df+x^2)^((df+1)/2))
    lfac <- (df/2)*log(df) -.5*ncp^2 - (log(pi)/2 + lgamma(df/2)  +(df+1)/2* log(df+x^2))
    if(check) stopifnot(all.equal(fac, exp(lfac), tol=tol.check))
    k <- 0:M
    if(check) suppressWarnings(
    terms  <-  gamma((df + k + df)/2) * ncp^k    /  factorial(k) * (2*x^2/(df+x^2)) ^ (k/2)  )
    lterms <- lgamma((df + k + df)/2)+k*log(ncp) - lfactorial(k) + (k/2) * log(2*x^2/(df+x^2))
    if(check) {
        ii <- which(!is.na(terms))
        stopifnot(all.equal(exp(lterms[ii]), terms[ii], tol=tol.check))
    }
    lf <- lfac + lsum(lterms)
    if(log) lf else exp(lf)
}
dntRwrong <- Vectorize(dntRwrong1, c("x", "df", "ncp"))
}##-- just for historical reference

### Johnson, Kotz and Balakrishnan (1995) [2nd ed.] have
###  (31.15) [p.516] and (31.15'), p.519 -- and they  contradict by a factor  (1 / j!)
###  in the sum
### ==> via the following: "prove" that (31.15) is correct, and  (31.15') is missing the "j!"
.dntJKBch1 <- function(x, df, ncp, log = FALSE, M = 1000, check=FALSE, tol.check = 1e-7)
{
    stopifnot(length(x) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              ncp >= 0, df > 0,
              is.numeric(M), M >= 1, M == round(M))
    if(check) ## this is the formula (31.15) unchanged:
     fac <- exp(-.5*ncp^2) * gamma((df+1)/2) / (  sqrt(pi*df)* gamma(df/2)) *    (df/(df+x^2))^((df+1)/2)
    lfac <-     -.5*ncp^2  +lgamma((df+1)/2) - (.5*log(pi*df)+lgamma(df/2)) + log(df/(df+x^2))*((df+1)/2)
    if(check) stopifnot(all.equal(fac, exp(lfac), tol=tol.check))
    j <- 0:M
    if(check)  suppressWarnings(## this is the formula (31.15) unchanged:
        ## [note that 1/gamma((df+1)/2) is indep. of j!]
    terms  <-  gamma((df + j + 1)/2) / ( factorial(j) * gamma((df + 1)/2))*    (x*ncp*sqrt(2)/sqrt(df+x^2))^ j )
    lterms <- lgamma((df + j + 1)/2) - (lfactorial(j) +lgamma((df + 1)/2))+ log(x*ncp*sqrt(2)/sqrt(df+x^2))* j
    if(check) {
        ii <- which(!is.na(terms))
        stopifnot(all.equal(exp(lterms[ii]), terms[ii], tol=tol.check))
    }
    lf <- lfac + lsum(lterms)
    if(log) lf else exp(lf)
}
.dntJKBch <- Vectorize(.dntJKBch1, c("x", "df", "ncp"))

## Q{MM}: is the [ a < x ] cutoff really exactly optimal?
logr <- function(x, a) { ## == log(x / (x + a)) -- but numerically smart; x >= 0, a > -x
    if(length(aS <- a < x) == 1L) {
        if(aS) -log1p(a/x) else log(x / (x + a))
    } else { # "vectors" : do ifelse(aS, .., ..) efficiently:
        r <- a+x # of correct (recycled) length and type (numeric, mpfr,  ..)
        r[ aS] <- -log1p((a/x)[aS])
        r[!aS] <- log((x / r)[!aS])
        r
    }
}

## New "optimized" and  "mpfr-aware" and *vectorized* (!) version:
dntJKBf <- function(x, df, ncp, log = FALSE, M = 1000)
{
    stopifnot(length(M) == 1, df >= 0, is.numeric(M), M >= 1, M == round(M))
    ln2 <- log(2)
    ._1.1..M <- c(1L, seq_len(M)) # cumprod(.) = (0!, 1!, 2! ..) =  (1, 1, 2, 6, ...)
    isN <- is.numeric(x) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(x, df, ncp)
    if(!isN) {
        if(isMpfr) {
	    stopifnot(requireNamespace("Rmpfr"))
            mpfr <- Rmpfr::mpfr ; getPrec <- Rmpfr::getPrec
	    prec <- max(getPrec(x), getPrec(df), getPrec(ncp))
	    pi <- Rmpfr::Const("pi", prec = max(64, prec))
	    ## this *is* necessary / improving results!
	    if(!inherits(x,  "mpfr")) x   <- mpfr(x,  prec)
	    if(!inherits(df, "mpfr")) df  <- mpfr(df, prec)
	    if(!inherits(ncp,"mpfr")) ncp <- mpfr(ncp,prec)
	    ln2 <- log(mpfr(2, prec))
	    ._1.1..M <- mpfr(._1.1..M, prec)
        } else {
            warning(" Not 'numeric' but also not  'mpfr' -- untested, beware!!")
        }
        ## needed for printing mpfr numbers {-> pkg Rmpfr}, e.g.
	## .N <- if(isMpfr) Rmpfr::asNumeric else as.numeric
    }

    x2 <- x^2
    lfac <- -ncp^2/2  - (.5*log(pi*df)+lgamma(df/2)) + logr(df, x2)*(df+1)/2
    j <- 0:M
    lfact.j <- cumsum(log(._1.1..M)) ## == lfactorial(j)
    nd <- length(df)
    LogRt <- 2*log(abs(ncp)) + ln2 + logr(x2, df) # (full length)
    ## now vectorize "lSum(x,df,ncp)" :
    lSum <- dx <- ncp*x + 0*df # delta * x  [of full length],  (correct == 0 for dx == 0)
    for(i in seq_along(dx)) { # --- compute lSum[i] ----------------
        alt <- dx[i] < 0 ## if(alt)  alternating sum !
	## use abs(ncp) : if(ncp < 0) ncp <- -ncp
	##lterms <- lgamma((df+j + 1)/2) - lfact.j + log(x*ncp*sqrt(2)/sqrt(df+x^2))* j
	lterms <- lgamma((df[1L+ (i-1L)%% nd] + j + 1)/2) - lfact.j + LogRt[i] * j/2
        lSum[i] <-
            if(alt) ## this is hard: even have *negative* sum {before log(.)} with mpfr,
                ## e.g. in  dnt.1(mpfr(-4, 128), 5, 10)   ???
                lssum(lterms, signs = c(1,-1), strict=FALSE)
            else
                lsum(lterms)
    }
    lf <- lfac + lSum
    if(log) lf else exp(lf)
}
## No longer, as have vectorized above!
## dntJKBf <- Vectorize(dntJKBf1, c("x", "df", "ncp"))
## instead, from 2019-10-04 :
dntJKBf1 <- function(x, df, ncp, log = FALSE, M = 1000) {
    .Deprecated("dntJKBf")
    dntJKBf(x=x, df=df, ncp=ncp, log=log, M=M)
}


## Orig: ~/R/MM/NUMERICS/dpq-functions/noncentral-t-density-approx_WV.R
##
## From: Wolfgang Viechtbauer <wviechtb@s.psych.uiuc.edu>
## To: <r-help@stat.math.ethz.ch>
## Subject: Re: [R] Non-central distributions
## Date: Fri, 18 Oct 2002 11:09:57 -0500 (CDT)
## .......
## This is an approximation based on Resnikoff & Lieberman (1957).
## .. quite accurate. .....
## .......
## .......
##
## MM: added 'log' argument and implemented  log=TRUE
## --- TODO: almost untested by MM [but see Wolfgang's notes in *_WV.R (s.above)]
dtWV <- function(x, df, ncp=0, log=FALSE) {
   dfx2 <- df + x^2 # = 'f+t^2' in Resnikoff+L.(1957), p.1 (by MM)
   y <- -ncp*x/sqrt(dfx2) # = 'y' in R.+L., p.1
   ## MM(FIXME): cancellation for y >> df  here :
   a <- (-y + sqrt(y^2 + 4*df)) / 2 # NB a = 't' in R.+L., p.25
   dfa2 <- df+a^2 ## << MM(2)
   if(log) {
       lHhmy <- df*log(a) + -0.5*(a+y)^2 +
           0.5*log(2*pi*a^2/dfa2) +
           log1p( - 3*df/(4*dfa2^2) + 5*df^2/(6*dfa2^3))
       lHhmy - (((df-1)/2)*log(2) + lgamma(df/2) + .5*log(pi*df)) +
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


###-- qnt() did not exist yet at the time I wrote this ...
##    ---
qtAppr <- function(p, df, ncp, lower.tail = TRUE, log.p = FALSE,
                    method = c("a","b","c"))
{
  ## Purpose: Quantiles of approximate non-central t
  ##  using Johnson,Kotz,.. p.521, formula (31.26 a) (31.26 b) & (31.26 c)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?qt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date:  6 Feb 99
    method <- match.arg(method)

  ##----------- NEED df >> 1 (this is from experiments below; what exactly??)
  z <- qnorm(p, lower.tail=lower.tail, log.p=log.p)
  if(method %in% c("a","c")) {
      b <- b_chi(df)
      b2 <- b*b
  }
  ## For huge `df';  b2 ~= 1 ---> method b) sets b = b2 = 1
  switch(method,
         "a" = {
             den <- b2 - z*z*(1-b2)
             (ncp*b + z*sqrt(den + ncp^2*(1-b2)))/den
         },
         "b" = {
             den <- 1 - z*z/(2*df)
             (ncp + z*sqrt(den + ncp^2/(2*df)))/den

         },
         "c" = {
             den <- b2 - z*z/(2*df)
             (ncp*b + z*sqrt(den + ncp^2/(2*df)))/den
         })
}
