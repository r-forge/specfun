#### R simulation of C's pnchisq_raw()  -*- delete-old-versions: never -*-
#### ---------------------------------

## For using this, see ./chisq-nonc-ex.R , for testing  ./pnchisq-tst.R
##		       ~~~~~~~~~~~~~~~~~                ~~~~~~~~~~~~~~~

### Exports:
### -------
### pnchisq  (q, df, ncp = 0, errmax = 1e-12, itrmax = 10* 10000, verbose = 1)
### rr       (i, lambda)
### titleR.exp  [expression]
### plotRR   (lambda, iset = 1:(2*lambda), do.main=TRUE)

## ~/R/D/r-patched/R/src/nmath/dpq.h :
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")
## for the R.D*() functions and vars :

pnchisq <- function(q, df, ncp = 0, lower.tail = TRUE,
                    ## log.p = FALSE,
                    cutOffncp = 80, itSimple = 110,
                    errmax = 1e-12, reltol = 1e-11,
                    itrmax = 10* 10000,
                    verbose = 1, xLrg.sigma = 5)
{
  ## Purpose: R simulation of C's pnchisq_raw()
  ## ----------------------------------------------------------------------
  ## Arguments: as pchisq()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14 Apr 2003, 21:17

    ## Use this such that verbose <= 0 is 'quiet'
    Cat <- function(...) if(verbose > 0) cat(...)

    if(length(q) != 1 || length(df) != 1 || length(ncp) != 1)
        stop("arguments must have length 1 !")
    ## follow C code :
    x <- q
    if(x <= 0) {
	if(x == 0 && df == 0)
	    return(if(lower_tail) exp(-0.5*ncp) else -expm1(-0.5*ncp))
        ## else
        return(R.DT.0(lower.tail, log.p=FALSE))
    }
    if(!is.finite(x))
        return(R.DT.1(lower.tail, log.p=FALSE))
    lam <- 0.5 * ncp

    if(ncp < cutOffncp) {
        ## C code uses LDOUBLE here ..
        pr <- exp(-lam)
        Cat(sprintf("pnchisq(x=%10g, ncp < cutoff): pr = %g ..", x, pr))
	sum <- sum2 <- 0
	for(i in 0:(itSimple-1)) {
            sum2 <- sum2 + pr
	    sum <- sum+ pr * pchisq(x, df+2*i,
                                    lower.tail=lower.tail, log.p=FALSE)
	    if (sum2 >= 1-1e-15) break
            pr <- pr * lam/(i <- i+1)
        }
        Cat(sprintf(" ==> final pr=%g, i = %3d, sum2 = %18.16f\n", pr, i, sum2))
        return(sum/sum2)
    }
    ## else :

    dbl.min.exp <- log(2) * .Machine$double.min.exp
    ##/*= -708.3964 for IEEE double precision */
    lamSml <-  (-lam < dbl.min.exp)
    if (lamSml) {
        ##  stop("non centrality parameter (= ",lam,
        ##       ") too large for current algorithm")
        Cat("large 'lambda' ==> working with ln(u)\n")
        u <- 0
        lu <- -lam # == ln(u)
        l.lam <- log(lam)
    }
    else
        u <- exp(-lam)
    v <- u
    x2 <- 0.5* x
    f2 <- 0.5* df
    fx.2n <- df - x ##  f - x

    if(f2 * .Machine$double.eps > 0.125 &&
       abs(t <- x2 - f2) < sqrt(.Machine$double.eps) * f2) {
	##/* evade cancellation error */
	## t <-  exp((1 - t)*(2 - t/(f2 + 1))) / sqrt(2* pi*(f2 + 1))
        lt <- (1 - t)*(2 - t/(f2 + 1)) - 0.5 * log(2* pi*(f2 + 1))
	Cat(" (case I) ==>")
    }
    else {
 	##/* Usually (case 2): careful not to overflow .. : */
	lt <- f2*log(x2) - x2 - lgamma(f2 + 1)
        ## FIXME: if f2 = df/2 << 1  ----> use  lgamm1p(f2) !
    }
    Cat(" lt=", formatC(lt))
    tSml <- (lt < dbl.min.exp)
    if (tSml) {
        Cat(" => exp(lt) underflow protection")
        if(x > df + ncp +  xLrg.sigma * sqrt( 2*(df + 2*ncp))) {
            ## x > E[X] + xL.. * sigma(X)
             Cat(" and x > E(X) + ",formatC(xLrg.sigma),
                 "*sigma(X) : too large --> 1 \n")
             return(R.DT.1(lower.tail, log.p=FALSE)) ## better than 0 --- but definitely "FIXME"
        } ## else :
        l.x <- log(x)
        Cat(" ln(x)=",l.x,"\n")
        ans <- term <- t <- 0
        ##-   if(t <= 0)
        ##-       warning("too small lt: too small/large x (=",formatC(x),
        ##-               ") or large centrality parameter ",formatC(ncp),
        ##-               " for current algorithm;\n",
        ##-               "\t result is doubtful")
    }
    else {
        t <- exp(lt)
        Cat(", t=exp(lt)=", formatC(t),"\n")
        ans <- term <- v * t
    }

    ##/* check if (f+2n) is greater than x */

    n <- 1
    f2n <- df + 2#/* = f + 2*n */
    fx.2n <- fx.2n + 2#;/* = f - x + 2*n */

    firstBound <- TRUE # only for debug-printing
## in C; rather use a for(n = 1, f2n = .., .;  <no test>; n++, f2n += 2, ...) {
    while(TRUE) {
	if(verbose >= 2) cat("_OL_: n=",n,"")

	## f_2n    === f + 2*n
	## f_x_2n  === f - x + 2*n   > 0  <==> (f+2n)  >   x
	if (fx.2n > 0) {
	    ##/* find the error bound and check for convergence */

            ## print only the first time! ==> need xtra variable firstBound
            if(firstBound) {
                Cat(" n= ",n,", fx.2n = ",formatC(fx.2n)," > 0\n",sep='')
                firstBound <- FALSE
            }
	    bound <- t * x / fx.2n
            ##ifdef DEBUG_pnch
            ## REprintf("\n L10: n=%d; term= %g; bound= %g",n,term,bound);
            is.r <- is.it <- FALSE
	    if (((is.b <- bound <= errmax) &&
                 (is.r <- term  <= reltol * ans)) || (is.it <- n > itrmax))
            {
                Cat("BREAK n=",n, if(is.it) "> itrmax",
                    "; bound= ",formatC(bound), if(is.b)"<= errmax",
                    "rel.err= ",formatC(term/ans),if(is.r)"<= reltol\n")
                break                   # out completely
            }

	}

        ##/* evaluate the next term of the */
        ##/* expansion and then the partial sum */

        if(lamSml) {
            lu <- lu + l.lam - log(n) ## u <- u* lam / n
            if(lu >= dbl.min.exp) {
                ## no underflow anymore ==> change regime
                Cat("  n=",n,
                    "; nomore underflow in u = exp(lu) ==> change\n")
                v <- u <- exp(lu) ## the first non-0 'u'
                lamSml <- FALSE
            }
        } else {
            u <- u* lam / n
            v <- v + u
        }
        if(tSml) {
            lt <- lt + l.x - log(f2n) ## t <- t * (x / f2n)
            if(lt >= dbl.min.exp) {
                ## no underflow anymore ==> change regime
                Cat("  n=",n,
                    "; nomore underflow in t = exp(lt) ==> change\n")
                t <- exp(lt) ## the first non-0 't'
                tSml <- FALSE
            }
        } else {
            t <- t * (x / f2n)
        }
        if(!lamSml && !tSml) {
            term <-  v * t
            if(verbose >= 2)
                cat(" il: term=",formatC(term,wid=10),
                    "rel.term=",formatC(term/ans, wid=10),"\n")
            ans <- ans + term
        } else if(verbose >= 2) cat(".")
        n <- n+1
        f2n <- f2n + 2
        fx.2n <- fx.2n + 2

    } ## while(TRUE)

    ## L_End:
    if (bound > errmax) { ## NOT converged
	warning("pnchisq(x,....): not converged in ",itrmax," iter.")
    }
    ##ifdef DEBUG_pnch
    ## REprintf("\n == L_End: n=%d; term= %g; bound=%g\n",n,term,bound);
    structure(R.D.Lval(ans, lower.tail=lower.tail), iter = n)
}

## Cheaply Vectorized version:
pnchisqV <- function(x, ..., verbose = 0)
    sapply(x, pnchisq, ..., verbose = verbose)

pnchisq.Pat <- function(q, df, ncp = 0,
                       lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Patnaik(1949)'s approximation to pnchisq() - using pchisq()
    ##    This is also the one in Abramowitz & Stegun, p.942, 26.4.27
    ## Johnson,Kotz,...: This is  O(1/sqrt(ncp)) for ncp -> Inf
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Feb 2004, 15:32
    e <- df + ncp
    d <- e + ncp ## = df + 2*ncp
    ic <- e/d
    pchisq(q*ic, df=e*ic, lower.tail=lower.tail, log.p=log.p)
}

pnchisq.Pea <- function(q, df, ncp = 0,
                       lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Pearson(1959) approximation to pnchisq() - using pchisq()
    ## Johnson,Kotz,...: Error is O(1/ncp) for ncp -> Inf --
    ## 		i.e. better than Patnaik for right tail
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Feb 2004, 15:38
    n2 <- df + 2*ncp
    n3 <- n2 + ncp ## = df + 3*ncp
    r <- n2 / n3
    pchisq((q + ncp*ncp/n3) *r, df = n2*r*r,
           lower.tail=lower.tail, log.p=log.p)
}

## Abdel-Aty (1954) .. Biometrika
pnchisq.AbdelAty <-  function(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Abdel-Aty(1954) "first approx." (Wilson-Hilferty) approximation to pnchisq()
    ## Johnson,Kotz,...: Has it *WRONGLY*  in eq. (29.61a), p.463
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 2018-08-16
    r <- df + ncp
    n2 <- r + ncp # = df + 2*ncp
    V <- 2*n2 / (9*r^2)
    pnorm((q/r)^(1/3), mean = 1-V, sd = sqrt(V),
          lower.tail=lower.tail, log.p=log.p)
}

pnchisq.Sankaran.d <- function(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Sankaran(1959,1963) according to Johnson et al. p.463, (29.61d):
    r <- df + ncp  #              = ν + λ
    n2 <- r + ncp  # = df + 2*ncp = ν + 2λ
    n3 <- n2 + ncp # = df + 3*ncp = ν + 3λ
    q1 <- n2/r^2   # = (n + 2λ)/(n+λ)²
    h1 <- - 2/3 * r*n3/n2^2 # =  h - 1
    h <- 1 + h1
    mu <- 1 + h*h1*q1*(1 + (h-2)*(1-3*h)*q1/2)
    V  <-   2*h^2 *q1*(1 +  h1 * (1-3*h)*q1)
    pnorm((q/r)^h, mean = mu, sd = sqrt(V),
          lower.tail=lower.tail, log.p=log.p)
}

rr <- function(i, lambda)
{
    ## Purpose: r_v(i) :=  (v^i / i!) / e_{i-1}(v), where
    ##      e_n(x) := 1 + x + x^2/2! + .... + x^n/n! (n-th partial of exp(x))
    ## As function of i :
    ## 	o  Can this be put in a simple formula?
    ##  o  When is it maximal?
    ##  o  When does r_{lambda}(i) become smaller than (f+2i-x)/x = a + b*i ?
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Jan 2004, 18:13
    if(lambda < -log(2)*.Machine$double.min.exp)
        dpois(i, lambda) / ppois(i-1, lambda)# << does vectorize in i
    else # large lambda
        exp(i*log(lambda) - lgamma(i+1) +
            -lambda- ppois(i-1, lambda, log = TRUE))
}

titleR.exp <-
    expression(rr[lambda](i) ==
        frac(lambda^i / i*"!",
             1+ lambda+ lambda^2/2*"!" +cdots+ lambda^{i-1}/(i-1)*"!"))

plotRR <- function(lambda, iset = 1:(2*lambda), do.main=TRUE,
                   log = 'xy', cex = 0.4, col = c("red","blue"))
{
    ii <- sort(iset)
    if(do.main) {
        pmar <- par("mar"); on.exit(par(mar=pmar))
        par(mar=pmar + c(0,0,1.2,0))
    }
    plot(ii, rr(ii, lambda=lambda), log = log, cex = cex, col = col[1],
         type = 'o', xlab = "i", ylab = "r(i)",
         sub = substitute(lambda==l, list(l=lambda)),
         main = if(do.main) titleR.exp)
    lines(ii, lambda/ii, col = col[2])
    legend(ii[length(ii)], lambda, expression(rr[lambda](i), lambda / i),
           col = col, lty = 1, pch = c(1,NA), xjust = 1, bty = 'n')
}




###--- Consider the sums directly :
##
## Approach used in pnchisq.c for  ncp < 80 :
pnchisqTerms <-  function(x, df, ncp, lower.tail = TRUE, i.max = 1000)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  1 Apr 2008, 17:52

    if(length(x)   > 1) stop("'x' must be of length 1")
    if(length(df)  > 1) stop("'df' must be of length 1")
    if(length(ncp) > 1) stop("'ncp' must be of length 1")
    stopifnot(length(i.max) == 1, i.max > 1)

    lambda <- 0.5 * ncp
    sum <- 0.
    pr <- exp(-lambda)
    k <- 1:i.max
    p.k <- dpois(k, lambda) ## == pr * cumprod(lambda / k)
    t.k <- pchisq(x, df + 2*(k-1), lower.tail=lower.tail, log.p=FALSE)
    s.k <- p.k * t.k
    kMax <- which.max(s.k)
    list(p.k=p.k, t.k=t.k, s.k = s.k, kMax = kMax,
         f = sum(s.k[1:kMax]) + sum(s.k[i.max:(kMax+1L)]))
}## {pnchisqTerms}


## Approach for ...
ss <- function(x, df, ncp, i.max = 10000)
{
    ## Purpose:
    ## pnchisq() :=                         sum_{i=0}^{n*} v_i  t_i
    ##            = exp(-lambda) t_0(x,f) * sum_{i=0}^{n*} v'_i t'_i
    ## hence the summands
    ## s_i := v'_i * t'_i
    ## 		     t'_0 := 1 ;
    ##               t'_i = t'_{i-1} * x /(f+2i)
    ##                    = x^i / prod_{j=1}^i (f + 2j)
    ##
    ##        v'_0 := u'_0 = 1;
    ##        v'_k := v'_{k-1} + u'_k
    ##                           u'_k := u'_{k-1} * (lambda / i)
    ## Return: list(s =(s_i){i= 0:i2},  i1, max)
    ##		where i1  : the index of the first non-0 entry in s[]
    ##                i2  : in 0:i.max such that trailing 0's are left away
    ##                max : the (first) index of the maximal entry s[]
    ## ----------------------------------------------------------------------
    ## Arguments: as p[n]chisq(), but only scalar
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  6 Feb 2004, 19:27
    if(length(x)   > 1) stop("'x' must be of length 1")
    if(length(df)  > 1) stop("'df' must be of length 1")
    if(length(ncp) > 1) stop("'ncp' must be of length 1")
    expMiMa <- log(2)*unlist(.Machine[c("double.min.exp","double.max.exp")])
    ## ~= (-708.4, 709.8)  (for IEEE FP)

    lambda <- ncp/2

    i <- 1:i.max #- really: from 0 to i.max

    xq <- x / (df + 2*i)
    tt <- cumprod(c(1, xq))
    it0 <- tt == 0 # underflow
    useLt <- any(it0) ## useLt := use log terms in formula
    if(useLt) { ## work with log(tt) there
        ltt <- cumsum(c(0, log(xq)))
    }

    ## Nota Bene: v[n] == e_n(lambda) * exp(-lambda)  is always in [0,1)
    useLv <- !(expMiMa[1] < -lambda && 1/lambda < expMiMa[2])
    if(!useLv) {## otherwise overflows/underflows
        u <- exp(-lambda)*cumprod(c(1, lambda / i))
        v <- cumsum(u)
    } else { ## lambda quite large or small --> compute log(u)
        ##lu <- cumsum(c(-lambda, log(lambda / i)))
        ##lu <- cumsum(c(-lambda, log(lambda) - log(i)))
        ## lu + lambda:
        luPl <- c(0,i*log(lambda)) + cumsum(c(0, - log(i)))
        useLv <- useLt || any(luPl - lambda < expMin)
        if(useLv) {
            lv <- -lambda + log(cumsum(exp(luPl)))
            if(!useLt)
                ltt <- cumsum(c(0, log(xq)))
        }
        else
            v <- cumsum(exp(luPl -lambda))
    }
    ## the sequence  r := tt * v  :
    if(useLv)
        r <- exp(ltt + lv)
    else {
        if(useLt && any(it0 <- it0 & v != 0))
            tt[it0] <- exp(ltt[it0])
        r <- tt * v
    }

    ## now get it's attributes:

    d <- diff(p <- r > 0)
    i1 <- which.max(d) # [i1] -> [i1+1]: first change from 0 to >0
    i2 <- i.max+1 - which.min(rev(d))
    ## [i2] -> [i2+1]: last change from >0 to 0
    r <- r[1:i2]
    return(list(s = r, i1 = i1, max = which.max(r)))
}

ss2 <- function(x, df, ncp, i.max = 10000, eps = .Machine$double.eps)
{
    ## Purpose: "Statistic" on ss() ==> give only interesting indices
    ## Author: Martin Maechler, Date:  7 Feb 2004.
    sss <- ss(x,df,ncp,i.max)
    s <- sss$s
    i.need <- range(which(s > eps * s[sss$max]))
    c(i1=sss$i1, i2=length(s), iN = i.need, max = sss$max)
}


pnchisq.ss <- function(x, df, ncp, i.max = 10000) {
    ## Using ss() for the non-central chisq prob.
    si <- ss(x=x, df=df, ncp=ncp, i.max = i.max)
    ## old version: had exp(-ncp/2)*
    2*dchisq(x, df = df +2) * sum(si$s)
}

## Instead of limited (overflow!) ss(),
## use C - code which parallels  pnchisq()'s in C:
dyn.load("/u/maechler/R/MM/NUMERICS/dpq-functions/pnchisq-it.so")

pnchisqIT <- function(q, df, ncp = 0, errmax = 1e-12,
                      reltol = 2*.Machine$double.eps, itrmax = 1e5)
{
    if(length(q) != 1 || length(df) != 1 || length(ncp) != 1)
        stop("arguments must have length 1 !")
    r <- .C("Pnchisq_it",
            x = as.double(q),
            f = as.double(df),
            theta = as.double(ncp),
            errmax = as.double(errmax),
            reltol = as.double(reltol),
            itrmax = as.integer(itrmax),
            i0 = integer(1),
            n.terms = integer(1),
            terms = double(itrmax+1), ## !!
            prob  = double(1),
            DUP = FALSE)
    length(r$terms) <- r$n.terms
    r[c("prob", "i0", "n.terms", "terms")]
}

ss2. <- function(q, df, ncp = 0, errmax = 1e-12,
                 reltol = 2*.Machine$double.eps, itrmax = 1e5,
                 eps = reltol)
{
    ## Purpose: "Statistic" on ss() ==> give only interesting indices
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 16 Feb 2004, 10:57

    if(length(q) != 1 || length(df) != 1 || length(ncp) != 1)
        stop("arguments must have length 1 !")
    r <- .C("Pnchisq_it",
            x = as.double(q),
            f = as.double(df),
            theta = as.double(ncp),
            errmax = as.double(errmax),
            reltol = as.double(reltol),
            itrmax = as.integer(itrmax),
            i0 = integer(1),
            n.terms = integer(1),
            terms = double(itrmax+1), ## !!
            prob  = double(1),
            DUP = FALSE)
    length(r$terms) <- r$n.terms
    nT <- length(s <- r$terms)

    ## now get it's attributes:
    if(nT > 1) {
        d <- diff(p <- s > 0)
        if(sum(diff(sign(d))) > 1) warning("more than local extremum")
        i1 <- which.max(d)     # [i1] -> [i1+1]: first change from 0 to >0
        i2 <- nT+1 - which.min(rev(d))## [i2] -> [i2+1]:last change from >0 to 0
        ## expect i2 == n.terms == nT :
        if(i2 != nT) { warning("i2 == ",i2," != nT = ",nT) ; s <- s[1:i2] }

        iMax <- which.max(s)
        i.need <- range(which(s > eps * s[iMax]))

        if(i2 != i.need[2] ## << happens:  then i2 == iN[2] + 1:
           &&
           i2 != i.need[2]+1) warning("i2 == ",i2," != iN[2] = ",i.need[2])

    } else { ## only 1 term
        i1 <- i2 <- NA
        iMax <- 1
        i.need <- c(1,1)
    }
    c(i0 = r$i0, nT = r$n.terms, i1=i1, i2=i2, iN = i.need, iMax = iMax)
}
