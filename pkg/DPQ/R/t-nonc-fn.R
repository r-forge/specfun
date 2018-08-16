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
###        -----  is very similar to  b.nu (and can use same asymptotic)

b.nu <- function(nu, one.minus = FALSE, c1 = 341, c2 = 1000)
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
    r[BL] <- b.nu.asymp(nu[BL], one.minus=one.minus)
    ##       ----------
    r
}

b.nu.asymp <- function(nu, one.minus = FALSE)
{
  ## Purpose: Compute  E[ chi_nu ] / sqrt(nu)  --- for "LARGE" nu
  ## ----------------------------------------------------------------------
  ## Arguments: nu >=0  (degrees of freedom)
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 11 Mar 1999, 15:11

  ## Abramowitz & Stegun, 6.1.47 (p.257) for a= 1/2, b=0 : __ for  z --> oo ___
  ## gamma(z + 1/2) / gamma(z) ~ sqrt(z)*(1 - 1/(8z) + 1/(128 z^2) + O(1/z^3))
  ## i.e. for b(nu), z = nu/2;  b(nu) = sqrt(1/z)* gamma((nu+1)/2) / gamma(nu/2)
  ##  b(nu) = 1 - 1/(8z) + 1/(128 z^2) + O(1/z^3) ~ 1 - 1/8z * (1 - 1/(16z))
  qq <- 1/(8*nu) # = 1/(16z)
  qq <- 2*qq*(1- qq)
  if(one.minus) qq else 1-qq
}


pt.appr1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE) {
    ## Approximate pt()  [non-central t] --- using the "extreme case" formula
    ##  from C's pnt.c  and  pntR() below:
    ##
    ## if (df > 4e5 || del*del > 2*M_LN2*(-(DBL_MIN_EXP))) {
    ## /*-- 2nd part: if del > sqrt(2*log(2)*1021) = 37.62.., then p=0 below
    ## FIXME: test should depend on `df', `tt' AND `del' ! */
    ## /* Approx. from	 Abramowitz & Stegun 26.7.10 (p.949) */

    ## VECTORIZE this:
    ## if(t >= 0) {
    ##     negdel <- FALSE
    ##     del <- ncp
    ## } else {
    ##     negdel <- TRUE
    ##     t   <- -t
    ##     del <- -ncp
    ## }
    negdel <- t < 0
    t  [negdel] <- -   t[negdel]
    del[negdel] <- - del[negdel]

    s <- 1/(4*df)
    pnorm(t*(1 - s), m = del, s = sqrt(1 + t*t*2*s),
          lower.tail = (lower.tail != negdel), log.p=log.p)
}

pt.appr <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE)
{
  ## Purpose: approximate non-central t -- still works FAST for huge ncp
  ##      [but has *wrong* asymptotic tail (for |t| -> oo)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?pt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date:  6 Feb 1999, 11 Mar 1999
  ## Johnson, Kotz,... , p.520, after (31.26a)

  b <- b.nu(df)
  ##   =====
  pnorm((t*b - ncp)/sqrt(1+ t*t*(1 - b*b)),
        lower.tail = lower.tail, log.p = log.p)
}

c.nu <- function(nu) {
    ## Purpose: The constant in  dt(t, nu, log=TRUE) =
    ##         = \log f_{\nu}(t) = c_{\nu} - \frac{\nu+1}{2}\log(1 + x^2/\nu)
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  1 Nov 2008, 14:26

    ## Limit nu -> Inf: c_{\nu} \to - \log(2\pi)/2 = -0.9189385
    r <- nu
    nu. <- nu[I <- nu < 200]
    r[I] <- log(gamma((nu.+1)/2)/gamma(nu./2)/sqrt(pi*nu.))
    nu. <- nu[I <- !I & nu < 1e4]
    r[I] <- lgamma((nu.+1)/2) - lgamma(nu./2) - log(pi*nu.)/2
    I <- nu >= 1e4
    r[I] <- c.nu.asymp(nu[I])
    r
}

c.nu.asymp <- function(nu)
{
    ## Purpose: Asymptotic of c.nu -- via Abramowitz & Stegun, 6.1.47 (p.257)
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## -log(2*pi)/2 -1/(4*nu) * (1 + 5/(96*nu^2))
    ##                               ^^^^^^^^^^^ not quite ok
    -log(2*pi)/2 -1/(4*nu)
}

c.pt <- function(nu)
{
  ## Purpose: the asymptotic constant in log F_nu(-x) = const(nu) - nu * log(x)
  ##          where F_nu(x) == pt(x, nu)
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 Nov 2008, 09:34
  c.nu(nu) + (nu-1)/2 * log(nu)
}


## MM:  My conclusion of experiments in ./pnt-precision-problem.R :
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
              df > 0, ncp > 0) ## ncp == 0 --> pt()

    ## Make this workable also for "mpfr" objects --> use format() in cat()
    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    if(!isN) {
	if(verbose) cat("some 'non-numeric arguments .. fine\n")
	if((isMpfr <- is(t, "mpfr") || is(df, "mpfr") || is(ncp, "mpfr"))) {
	    stopifnot(require("Rmpfr"))
	    ## if(!exists("pbetaRv1", mode="function"))
	    ##     source("~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/pbetaR.R")
	    prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
	    pi <- Const("pi", prec = max(64, prec))
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
        pnt.. <<- pnorm(tt*(1 - s), del, sqrt(1. + tt*tt*2.*s),
                        lower.tail = (lower.tail != negdel), log.p=log.p)
        cat("large 'df' or \"large\" 'ncp' -- C code would return pnorm(*) =: pnt.. =",
            format(pnt.., digits=16), "\n")
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
                    if(verbose == 1) sprintf("%d iter.: ", it) else "\\--> ",
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
##' @author Marius Hofert, Martin Maechler
lsum <- function(lx, l.off = max(lx)) l.off + log(sum(exp(lx - l.off)))
## FIXME: the case with lx == -Inf for all should give -Inf, but gives NaN
## ==> fixed in copula pkg (R-forge)


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
##' @author Marius Hofert and Martin Maechler
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
    if(!isN) {
        if((isMpfr <- is(t, "mpfr") || is(df, "mpfr") || is(ncp, "mpfr"))) {
            stopifnot(require("Rmpfr"))
            ## if(!exists("pbetaRv1", mode="function"))
            ##     source("~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/pbetaR.R")
            prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
            pi <- Const("pi", prec = max(64, prec))
            pbeta <- Vectorize(pbetaRv1, "qin") # so pbeta(x, p, <vector q>) works
            dbeta <- function(x, a,b, log=FALSE) {
                lval <- (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b)
                if(log) lval else exp(lval)
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
            R.Log1.Exp(-LS) ## = log1mexp(LS)
        else ## upper tail: log(1 - (1 - exp(-LS))) = -LS
            -LS
    } else {
        if(lower.tail) -expm1(-LS) else exp(-LS)
    }
}
pnt3150  <- Vectorize(pnt3150.1, c("t", "df", "ncp"))

### New version of pntR1(), pntR() ----  Using the  Posten (1994) algorithm
pntP94.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                   itrmax = 1000, errmax = 1e-12, verbose = TRUE) {

    stopifnot(length(t) == 1, length(df) == 1, length(ncp) == 1,
              df > 0, ncp > 0) ## ncp == 0 --> pt()

    Cat <- function(...) if(verbose > 0) cat(...)
    ## Make this workable also for "mpfr" objects --> use format() in cat()
    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    if(!isN) {
        if((isMpfr <- is(t, "mpfr") || is(df, "mpfr") || is(ncp, "mpfr"))) {
            stopifnot(require("Rmpfr"))
            ## if(!exists("pbetaRv1", mode="function"))
            ##     source("~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/pbetaR.R")
            prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
            pi <- Const("pi", prec = max(64, prec))
            pbeta <- pbetaRv1
            dbeta <- function(x, a,b, log=FALSE) {
                lval <- (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b)
                if(log) lval else exp(lval)
            }
        }
        ## needed for printing mpfr numbers {-> pkg Rmpfr}, e.g.
	.N <- if(isMpfr) asNumeric else as.numeric
        Cat("Some 'non-numeric arguments .. fine\n")
    }

    neg.t <- (t < 0)
    if(neg.t) {
        if(verbose) cat("t < 0  ==> swap delta := -ncp\n")
	tt <- -t; del <- -ncp
    } else {
        tt <- t; del <- ncp
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
            R.Log1.Exp(-LS) ## = log1mexp(LS)
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


## >> Just for historical reference :
##    =============================
dntR.1 <- function(x, df, ncp, M = 1000, log = FALSE, check=FALSE, tol.check = 1e-7)
{
    ## R's source ~/R/D/r-devel/R/src/nmath/dnt.c  claims -- from 2003 till 2014 -- but *WRONGLY*
    ## *	   f(x, df, ncp) =
    ## *		df^(df/2) * exp(-.5*ncp^2) /
    ## *		(sqrt(pi)*gamma(df/2)*(df+x^2)^((df+1)/2)) *
    ## *		sum_{k=0}^Inf  gamma((df + k + df)/2)*ncp^k /
    ## *				prod(1:k)*(2*x^2/(df+x^2))^(k/2)
    stopifnot(length(x) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              ncp >= 0, df >= 0,
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
dntR <- Vectorize(dntR.1, c("x", "df", "ncp"))

if(FALSE) {## no experiments in this file ... but for testing interactively :
dntR(1:6, df=3, ncp=5, check=TRUE)
## [1] 0.0032023 0.2728377 1.5174647 2.4481320 2.4106931 1.9481189 -- wrong "but sensible"

x <- seq(-1,20, by=1/4)
plot(x, dt(x, df=3, ncp=5) / dntR(x, df=3, ncp=5))
## clearly wrong

}## testing


### Johnson, Kotz and Balakrishnan (1995) [2nd ed.] have
###  (31.15) [p.516] and (31.15'), p.519 -- and they  contradict by a factor  (1 / j!)
###  in the sum
### ==> via the following: "prove" that (31.15) is correct, and  (31.15') is missing the "j!"
dnt.1 <- function(x, df, ncp, M = 1000, log = FALSE, check=FALSE, tol.check = 1e-7)
{
    stopifnot(length(x) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              ncp >= 0, df >= 0,
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
dnt <- Vectorize(dnt.1, c("x", "df", "ncp"))

if(FALSE) {## no experiments in this file ... but for testing interactively :
(ft <- dnt(1:6, df=3, ncp=5, check=TRUE))
## [1] 0.000508267 0.026060733 0.119137668 0.176104468 0.165771811 0.130541073
## correct !! *with* the factorial !
stopifnot(all.equal(ft, dt(1:6, df=3, ncp=5)))
}## testing

logr <- function(x, a) ## == log(x / (x + a)) -- but numerically smart; x > 0, a >= 0 > -x
    if(a < x) -log1p(a/x) else log(x / (x + a))

## New "optimized" and  "mpfr-aware" version:
dnt.1 <- function(x, df, ncp, M = 1000, log = FALSE, verbose=FALSE, tol.check = 1e-7)
{
    stopifnot(length(x) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              df >= 0, is.numeric(M), M >= 1, M == round(M))
    isN <- is.numeric(x) && is.numeric(df) && is.numeric(ncp)
    ln2 <- log(2)
    ._1.1..M <- c(1L, seq_len(M)) # cumprod(.) = (0!, 1!, 2! ..) =  (1, 1, 2, 6, ...)
    if(!isN) {
        if((isMpfr <- is(x, "mpfr") || is(df, "mpfr") || is(ncp, "mpfr"))) {
	    stopifnot(require("Rmpfr"))
	    prec <- max(getPrec(x), getPrec(df), getPrec(ncp))
	    pi <- Const("pi", prec = max(64, prec))
	    ## this *is* necessary / improving results!
	    if(!is(x,  "mpfr")) x   <- mpfr(x,  prec)
	    if(!is(df, "mpfr")) df  <- mpfr(df, prec)
	    if(!is(ncp,"mpfr")) ncp <- mpfr(ncp,prec)
	    ln2 <- log(mpfr(2,prec))
	    ._1.1..M <- mpfr(._1.1..M, prec)
        } else {
            warning(" Not 'numeric' but also not  'mpfr' -- untested, beware!!")
        }
        ## needed for printing mpfr numbers {-> pkg Rmpfr}, e.g.
	.N <- if(isMpfr) asNumeric else as.numeric
    }

    x2 <- x^2
    lfac <- -ncp^2/2  - (.5*log(pi*df)+lgamma(df/2)) + logr(df, x2)*(df+1)/2
    j <- 0:M
    lfact.j <- cumsum(log(._1.1..M)) ## == lfactorial(j)
    dx <- ncp*x # delta * x
    lSum <- if(dx == 0) 0 else {
	alt <- (sign(dx) == -1) ## if(alt)  alternating sum !
	if(ncp < 0) ncp <- -ncp
	##lterms <- lgamma((df+j + 1)/2) - lfact.j + log(x*ncp*sqrt(2)/sqrt(df+x^2))* j
	lterms <- lgamma((df + j + 1)/2) - lfact.j + (2*log(ncp) + ln2 + logr(x2, df)) * j/2
	if(alt) ## this is hard: even have *negative* sum {before log(.)} with mpfr,
            ## e.g. in  dnt.1(mpfr(-4, 128), 5, 10)   ???
            lssum(lterms, signs = rep_len(c(1,-1), length(j)), strict=FALSE)
        else
            lsum(lterms)
    }
    lf <- lfac + lSum
    if(log) lf else exp(lf)
}
dnt <- Vectorize(dnt.1, c("x", "df", "ncp"))


###-- qnt() did not exist yet at the time I wrote this ...
##    ---
qt.appr <- function(p, df, ncp, method = c("a","b","c"))
{
  ## Purpose: Quantiles of approximate non-central t
  ##  using Johnson,Kotz,.. p.521, formula (31.26 a) (31.26 b) & (31.26 c)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?qt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date:  6 Feb 99
    method <- match.arg(method)

  ##----------- NEED df >> 1 (this is from experiments below; what exactly??)
  z <- qnorm(p)
  if(method %in% c("a","c")) {
      b <- b.nu(df)
      b2 <- b*b
  }
  ## For huge `df';  b2 ~= 1 ---> method b) sets b = b2 = 1
  switch(method,
         "a" = {
             den <- b2 - z*z*(1-b2)
             (ncp*b + z*sqrt(den + ncp^2*(1-b2)))/den
         },
         "b" = {
             den <- 1 - z*z/(2*nu)
             (ncp + z*sqrt(den + ncp^2/(2*nu)))/den

         },
         "c" = {
             den <- b2 - z*z/(2*nu)
             (ncp*b + z*sqrt(den + ncp^2/(2*nu)))/den
         })
}
