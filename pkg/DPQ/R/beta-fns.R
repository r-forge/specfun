
##' Qab := Q(a,b) := Gamma(a+b) / Gamma(b)  =  Gamma(a) / Beta(a,b)
##' <==>
##' log(Q(a,b)) = logQab(a,b)
##'             = lgamma(a+b) - lgamma(b) =
##'             = lgamma(a)   - lbeta(a,b)
logQab_asy <- function(a,b, k.max = 5, give.all = FALSE)
{
  ## Purpose:  log(Q(a,b)) asymptotically when max(a,b) -> Inf  &  a^2/b -> C
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  8 Jun 97, 17:28

  if(any(a<=0) || any(b<=0)) stop("a & b must be > 0")
  if(all(a >= b)) { cc <- b; b <- a; a <- cc } #- Now: a <= b
  else if(any(a>b)) stop('Require a <= b  or  a > b   for  ALL vector elements')
  if(length(b)>1) stop("this function only works for SCALAR b (> a)")
  na <- length(a)
  pa <- rbind(Qab_terms(a, k.max)) ##  --> pa = matrix   na  x  (k.max+1)
  ##          ~~~~~~~~~
  ki <- if(k.max >= 1) k.max:1 else numeric(0)
  if(give.all) {
    r <- matrix(0, nrow=na, ncol=k.max+1,
                dimnames = list(format(a), paste(k.max:0)))
    for(k in ki) r[,k] <- pa[,k] + r[,k+1]*(a-k-1)/b
  } else {
    r <- numeric(na)
    for(k in ki) r <- pa[,k] + r*(a-k-1)/b
  }
  a*log(b) + log(1 + a*(a-1)/b * r)
}

Qab_terms <- function(a, k)
{
  ## Purpose: Compute the terms used for  Qab, and beta function
  ## 		Qab := Q(a,b) := Gamma(a+b) / Gamma(b)
  ## -------------------------------------------------------------------------
  ## Arguments: a: the smaller of the arguments of  beta(a,b)
  ##            k: the number of terms in the series expansion
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 19 Jun 97, 11:08

  k <- as.integer(k)
  if(!is.numeric(a)) a <- as.numeric(a) #- leave integers alone..
  if(length(k) != 1 || k > 5 || k < 0)
    stop("'k' must be 1 number in { 0, 1,..,5 }")
  na <- length(a)
  Qab.trms <- expression(
      rep.int(1/2, na),	##  k=1
      a/8 - 1/24, 	##  k=2
      a*(a-1)/48,	##  k=3
      (2/5 + a*(1+3*a*(a-2)))/1152,  ## k=4
      a*(1 + a*(2+a*(a-5)))/5760     ## k=5
  )[seq_len(k)]

  vapply(Qab.trms, eval, numeric(na), envir=environment())
}

### NOTA BENE: All these  lbeta*() are -- I think -- not better than the
##  ---------  builtin  lbeta()  [which exists only since 1998]

lbeta_asy <- function(a,b, k.max = 5, give.all = FALSE)
{
  ## Purpose:  log(beta( a, b))  which also works for HUGE  a / b or b/a
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ##  log(B(a,b)) = log(G(a) G(b) / G(a+b)) =
  ##              = log(G(a)) - log( G(a+b)/G(b) ) = log(G(a)) - log(Q(a,b))
  lgamma(a) - logQab_asy(a,b, give.all = give.all)
}


### FIXME ==> ../tests/qbeta-dist.R
##            -------- ~~~~~~~~~~~~ at the end has already *several* such approximations !!! <<< FIXME

##' Taylor Approximation to  log(a*B(a,b))  for  a << 1 -- derived from a --> 0
##'
##' for small a,  B(a,b) ~ 1/a, hence  a*B(a,b) ~ a * 1/a = 1, hence log(a*B(a,b)) ~= 0
##'
##' log(a * beta(a,b)) == log(a) + log(beta(a,b)) == log(a) + lbeta(a,b), but the latter suffers cancellation !!
##' We use the Taylor series
##' log(a * B(a,b)) = a \sum_{n=0}^\infty (\psi^{(n)}(1) - \psi^{(n)}(b)) \frac{a^n}{(n+1)!}
##'              where \psi^{(n)} is the n-th derivative of the "psi gamma" function
##' derived from Abramowitz & Stegun, p.259, (6.3.14) series for \psi(1+z)
##' it's integral is \log \Gamma(z+1) = - \gamma z + \sum_{k=2}^\infty \zeta(k)/k (-z)^k   (G)
##' and as \Gamma(z+1) = z \Gamma(z), have \log \Gamma(y) = \log(\Gamma(y+1)/y) = - \log(y) + formula (G) above.
##'
##' Now, \log(a \times B(a,b)) = \log(a) + \log \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}
##'                            = \log(a) + \log \Gamma(a)         + \log (\Gamma(b) / \Gamma(a+b))
##'                            = \gamma a + \pi^2/12 a^2 +O(a^3)  + \log (\Gamma(b) / \Gamma(a+b))
##'  ...
##'                            = a \times \sum_{n=0}^\infty .. sum above

##' @title Accurate Approximation of log(a * beta(a,b)) for small a
##' @param a, b  numeric vectors: the non-negative arguments of Beta
##' @return
##' @author Martin Maechler
laBeta <- function(a,b, nT) {
    stopifnot(length(nT) == 1, nT == as.integer(nT), nT >= 1)

### ===> ../tests/qbeta-dist.R
##       =========~~~~~~~~~~~~ at the end has already *several* such approximations !!! <<< FIXME

}


## MM: this is unused (why ??) -- should be "uniformly better" than lbeta_asy() above
lbetaMM <- function(a,b, cutAsy = 1e-2, verbose = FALSE)
{
  ## Purpose:  log(beta( a, b))  which also works for HUGE  a / b or b/a
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 May 97, 17:17
  if(a <= 0 || b <= 0) stop("a & b must be > 0")
  if(length(a)>1 || length(b) > 1)
    stop("this function only works for SCALAR arguments")
  stopifnot(is.numeric(cutAsy), length(cutAsy) == 1, cutAsy >= 0)
  if(a>b) { cc <- b; b <- a; a <- cc } #- Now: a <= b
  if(a*a < b*cutAsy) {
    if(verbose)
      message("a=",formatC(a)," b=",formatC(b), " -- using asymptotic  lbeta(.)")
    lgamma(a) - logQab_asy(a,b)
  } else
    lgamma(a)+lgamma(b)-lgamma(a+b) ## was =: lbeta00(a, b)
}

lbetaM <-  function(a,b, k.max = 5, give.all = FALSE)
{
  ## Purpose:  log(beta( a, b))  which also works for HUGE  a / b or b/a
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 May 97, 17:17
  if(any(a<=0) || any(b<=0)) stop("a & b must be > 0")
  na <- length(a)
  nb <- length(b)
  if(nb == 1 || nb == na) {
    r <- numeric(na)
    for(i in 1:na) r[i] <- lbeta_asy(a[i], b[if(nb>1)i else 1],
                                     k.max= k.max, give.all = give.all)
  } else { ##--- return  outer(..)
    r <- matrix(0, nrow = nb, ncol = na)
    for(i in 1:na)r[,i] <- lbeta_asy(a[i], b, k.max= k.max, give.all = give.all)
  }
  r
}
lbetaI <- function(a, n)
{
  ## Purpose:  log(beta(a, n)) for (small) INTEGER n.
  ##    beta(a,n) = G(a) G(n) / G(a+n) = (n-1)! / (a*(a+1)*...*(a+n-1))
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 May 97, 10:08
  if(is.na(n <- as.integer(n)) || min(n) <= 0 || max(n) > 10000)
    stop("'n' must be positive (not too large) integer")
  r <- numeric(nn <- length(n))
  for (i in 1:nn)
    r[i] <- lgamma(n[i]) - sum(log(a+0:(n[i]-1)))
  r
}

betaI <- function(a, n)
{
  ## Purpose: beta(a, n) for (small) INTEGER n.
  ##    beta(a,n) = G(a) G(n) / G(a+n) = 1*2*...*(n-1) / (a*(a+1)*...*(a+n-1))
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 May 97, 10:08
  if(length(n) > 1) stop("'n' must have length 1.")
  if(is.na(n <- as.integer(n)) || n <= 0)
    stop("'n' must be positive (not too large) integer")
  ni <- seq_len(n-1L)
  prod(ni) / prod(a+c(0,ni))
}

G_half <- 1.7724538509055160272981674833411451827974 # sqrt(pi) = Gamma(1/2)

if(FALSE) ## not closed form smart formula I think (2019-08)
lbetaIhalf <- function(a, n)
{
  ## Purpose:  log(beta(a, n + 1/2)) for (small) INTEGER n.
  ##    beta(a,n + h) = G(a) G(n+h) / G(a+n+h) = (n-1)! / (a*(a+1)*...*(a+n-1))
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 May 97, 10:08
  if(length(n) > 1) stop("'n' must have length 1.")
  if(is.na(n <- as.integer(n)) || n <= 0)
    stop("'n' must be positive (not too large) integer")

  stop("unfinished -- FIXME")
}

### pbeta()  --- is now used from TOM708 exclusively
### -------  ~/R/D/r-devel/R/src/nmath/toms708.c
###                                    ~~~~~~~~~

## From ~/R/D/r-devel/R/src/nmath/pgamma.c :
##                                --------
##/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */

## ---> lgamma1p() below

scalefactor <- 2^256 # == 2^2^2^3 == double prec. exact{checking via Rmpfr and gmp::as.bigz(2)^256}

## smallest exponent E for which exp(E) is "above underflow" {exp(E) == .Machine$double.xmin}
M_minExp <- log(2) * .Machine$double.min.exp # ~= -708.396

##/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ] =~=  -x */
M_cutoff <- log(2) * .Machine$double.max.exp / .Machine$double.eps ## = 3.196577e18


##' a useful 1-string  format(summary(<logical>)): __TODO: add to {sfsmisc} __
formatSummL <- function(x) {
    stopifnot(is.logical(x))
    if(!length(x)) return("logi(0)")
    ## else non-empty -- potentially has  {TRUE, FALSE, NA}
    tab <- table(x, exclude=NULL)
    if(any(I <- is.na(names(tab))))# "NA", not <NA>
        names(tab)[I] <- "NA"
    nt <- sort(names(tab), decreasing = TRUE)# typical starting with "TRUE"
    k <- length(
        snt <- if(any(x & !is.na(x))) sub("FALSE", "F.", nt, fixed=TRUE) else nt)
    m <- rbind(unname(tab[nt]), " ", snt
             , if(k >= 2) c(rep(", ", k-2), ", and ", "") # else NULL
             , deparse.level=0L)
    paste(m, collapse="")
}

## This is the "fast pure R" version that does iterations
## in parallel *where applicable*
## The next one, logcfR(), will entirely vectorize, ==> iterations separately
logcfR. <- function (x, i, d, eps, # ~ relative tolerance
                     maxit = 10000L, trace = FALSE)
{
    ## Continued fraction for calculation of  sum_{k=0}^Inf x^k/(i+k*d) =
    ##		1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ...
    ##
    ## auxiliary in log1pmx() and lgamma1p()
    stopifnot (i > 0, d >= 0, eps > 0,
               length(i) == 1, length(d) == 1, length(eps) == 1
               ## x < 1 for (i,d = 3,2) -- there x = -2 is fine, but in general?
               ## |x| < 1 # otherwise it does not (always?) converge,
               ##            ??
               )
    if(any_mpfr(x, i, d)) {
        stopifnot(requireNamespace("Rmpfr"))
        mpfr <- Rmpfr::mpfr ; getPrec <- Rmpfr::getPrec
        prec <- max(getPrec(x), getPrec(i), getPrec(d))
        mean <- Rmpfr::mean
        ## this *is* necessary / improving results!
        if(!inherits(x, "mpfr")) x <- mpfr(x, prec)
        if(!inherits(i, "mpfr")) i <- mpfr(i, prec)
        if(!inherits(d, "mpfr")) d <- mpfr(d, prec)
        asN <- Rmpfr::asNumeric
    } else asN <- as.numeric

    c1 <- 2 * d  # scalar
    c2 <- i + d  #  "
    c4 <- c2 + d #  "
    ## these will all depend on x :
    a1 <- c2
    b1 <- i * (c2 - i * x)
    r <- d * d * x
    A2 <- c4 * c2 - r
    B2 <- c4 * b1 - i * r

    it <- 0L
    if(trace) cat(sprintf("logcf(x[], i=%g, d=%g, eps=%g)", i, d, eps),
                  "  iterations:\n")
    n <- length(r <- x) # r[1:n] : the result (also good for 'mpfr' numbers)
    do.x <- rep(TRUE, n)
    while (any(needIt <- abs(A2 * b1 - a1 * B2) > abs(eps * b1 * B2))) {
        ## needIt, x, (a1,b1, A2, B2) are all of "shrinking length"
        ## (r, do.x) must stay at length == n
        it <- it+1L
        if(trace)
            cat(sprintf("it=%2d: needIt: %19s; ", it, formatSummL(needIt)))
        if(length(iFini <- which(!needIt))) {
            ## compute iFi, the indices in the *original* length-n vectors
            iFi <- which(do.x)[iFini]
            r   [iFi] <- A2[iFini] / B2[iFini]
            do.x[iFi] <- FALSE
            ## and shorten the rest:
            iK <- which(needIt) # the indices to keep
             x <-  x[iK]
            if(length(a1) == length(b1)) # (not when it = 1)
            a1 <- a1[iK]
            b1 <- b1[iK]
            A2 <- A2[iK]
            B2 <- B2[iK]
        }
        ## else : all(needIt) --> no finished parts, no shortening

        if(trace) cat(sprintf("length(x[<todo>])=%2d, ", sum(do.x)))

        c3 <- c2*c2 * x
	c2 <- c2 + d
	c4 <- c4 + d
	a1 <- c4 * A2 - c3 * a1
	b1 <- c4 * B2 - c3 * b1

	c3 <- c1*c1 * x
	c1 <- c1 + d
	c4 <- c4 + d
	A2 <- c4 * a1 - c3 * A2
	B2 <- c4 * b1 - c3 * B2

        m.B2 <- exp(mean(log(abs(B2))))
        if(trace) cat(sprintf(" m.B2=%12g", asN(m.B2)))
	if (m.B2 > scalefactor) {
            if(trace) cat(sprintf("  Lrg m.B2"))
	    a1 <- a1 / scalefactor
	    b1 <- b1 / scalefactor
	    A2 <- A2 / scalefactor
	    B2 <- B2 / scalefactor
	}
        else if (m.B2 < 1 / scalefactor) {
            if(trace) cat(sprintf("  Sml m.B2"))
	    a1 <- a1 * scalefactor
	    b1 <- b1 * scalefactor
	    A2 <- A2 * scalefactor
	    B2 <- B2 * scalefactor
	}
        if(trace) cat("\n")
        if(it > maxit) {
            warning("non-convergence in logcf(), iter > ", maxit)
            break
        }
    }
    if(trace && it <= maxit) cat("  logcf(*) end: after", it, "iterations.\n")

    r[do.x] <- A2 / B2
    ## and return
    r
} ## logcfR.()

## This is the entirely vectorizing "pure R" version
logcfR <- function (x, i, d, eps, # ~ relative tolerance
                    maxit = 10000L, trace = FALSE)
{
    ## Continued fraction for calculation of  sum_{k=0}^Inf x^k/(i+k*d) =
    ##		1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ...
    ##
    ## auxiliary in log1pmx() and lgamma1p()
    stopifnot (i > 0, d >= 0, eps > 0,
               length(i) == 1, length(d) == 1, length(eps) == 1
               ## x < 1 for (i,d = 3,2) -- there x = -2 is fine, but in general?
               ## |x| < 1 # otherwise it does not (always?) converge,
               ##            ??
               )
    if(any_mpfr(x, i, d)) {
        stopifnot(requireNamespace("Rmpfr"))
        mpfr <- Rmpfr::mpfr ; getPrec <- Rmpfr::getPrec
        prec <- max(getPrec(x), getPrec(i), getPrec(d))
        ## this *is* necessary / improving results!
        if(!inherits(x, "mpfr")) x <- mpfr(x, prec)
        if(!inherits(i, "mpfr")) i <- mpfr(i, prec)
        if(!inherits(d, "mpfr")) d <- mpfr(d, prec)
        asN <- Rmpfr::asNumeric
    } else asN <- as.numeric

    r <- x
    for(ii in seq_along(x)) {
        xi <- x[ii]
        c1 <- 2 * d
        c2 <- i + d
        c4 <- c2 + d
        a1 <- c2
        b1 <- i * (c2 - i * xi)
        ddx <- d * d * xi
        A2 <- c4 * c2 - ddx
        B2 <- c4 * b1 - i * ddx

        it <- 0L
        if(trace) cat(sprintf("logcf(x[%2d]=%8g, i=%g, d=%g, eps=%g)", ii, xi, i, d, eps),
                      if(trace >= 2)"  iterations:\n")
        while (abs(A2 * b1 - a1 * B2) > eps * abs(b1 * B2)) {
            c3 <- c2*c2*xi
            c2 <- c2 + d
            c4 <- c4 + d
            a1 <- c4 * A2 - c3 * a1
            b1 <- c4 * B2 - c3 * b1

            c3 <- c1 * c1 * xi
            c1 <- c1 + d
            c4 <- c4 + d
            A2 <- c4 * a1 - c3 * A2
            B2 <- c4 * b1 - c3 * B2

            it <- it+1L
            if(trace >= 2) cat(sprintf("it=%2d: ==> B2=%13g", it, asN(B2)))
            if (abs(B2) > scalefactor) {
                if(trace >= 2) cat(sprintf("  Lrg m.B2\n%27s", ""))
                a1 <- a1 / scalefactor
                b1 <- b1 / scalefactor
                A2 <- A2 / scalefactor
                B2 <- B2 / scalefactor
            }
            else if (abs(B2) < 1 / scalefactor) {
                if(trace >= 2) cat(sprintf("  Sml m.B2\n%27s", ""))
                a1 <- a1 * scalefactor
                b1 <- b1 * scalefactor
                A2 <- A2 * scalefactor
                B2 <- B2 * scalefactor
            }
            if(trace >= 2) cat(sprintf(" --> crit. |A2*b1 - a1*B2|/|b1*B2| = %g\n",
                                  asN(abs(A2 * b1 - a1 * B2)/abs(b1 * B2))))
            if(it > maxit) {
                warning("non-convergence in logcf(), iter > ", maxit)
                break
            }
        }
        if(trace && it <= maxit) cat(sprintf("  logcf(*) end: after %3d iterations.\n", it))
        r[ii] <- A2 / B2
        ##       -------
    } ## end for(ii ..)
    r
} ## logcfR() {"properly" vectorized}


logcf <- function (x, i, d, eps, trace = FALSE) {
    .Call(C_R_logcf, x, i, d, eps, trace)
}

## Accurate calculation of log(1+x)-x, particularly for small x.
## See also R-interface to R's C API [src/nmath/pgamma.c] log1pmx()  --> log1pmxC() in >> ./utils.R
## NB: minL1 was hard-wired 'minLog1Value'  in R's source of log1pmx()
log1pmx <- function(x, tol_logcf = 1e-14, eps2 = 0.01,
                    minL1 = -0.79149064, trace.lcf = FALSE,
                    logCF = if(is.numeric(x)) logcf else logcfR.)
{
    stopifnot(is.numeric(eps2), eps2 >= 0, is.numeric(minL1), -1 <= minL1, minL1 < 0)# < -1/4 ?
    r <- x
    if(any(c1 <- (x > 1 | x < minL1)))
        r[c1] <- log1p(x[c1]) - x[c1]
    ## else { ## ##/* expand in [x/(2+x)]^2 */
    if(any(c2 <- !c1)) {
        x <- x[c2]
	term <- x / (2 + x)
	y <- term * term
        r[c2] <- term * { ## not using ifelse(), rather what works with "mpfr"
            ## ifelse(abs(x) < eps2,
            ##        (((2 / 9 * y + 2 / 7) * y + 2 / 5) * y + 2 / 3) * y - x,
            ##        2 * y * logcf(y, 3, 2, tol_logcf) - x)
            A <- x
            if(any(isSml <- abs(x) <= eps2)) {
                i <- which(isSml)
                y. <- y[i]
                z <- if(is.numeric(x)) 2 else 2 + 0*y. # becomes y-like (e.g. mpfr)
                A[i] <- (((z / 9 * y. + z / 7) * y. + z / 5) * y. + z / 3) * y. - x[i]
            }
            if(length(iLrg <- which(!isSml))) {
                y. <- y[iLrg]
                A[iLrg] <- 2 * y. *
                    logCF(y., 3, 2, tol_logcf, trace = trace.lcf) - x[iLrg]
            }
            A
        }
    }
    r
}

lgamma1p <- function(a, tol_logcf = 1e-14, f.tol = 1., ...)
{
    ## Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5).

    eulers_const <- 0.5772156649015328606065120900824024

    ##/* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    N <- 40
    coeffs <-
        c(0.3224670334241132182362075833230126e-0, ##/* = (zeta(2)-1)/2 */
          0.6735230105319809513324605383715000e-1, ##/* = (zeta(3)-1)/3 */
          0.2058080842778454787900092413529198e-1,
          0.7385551028673985266273097291406834e-2,
          0.2890510330741523285752988298486755e-2,
          0.1192753911703260977113935692828109e-2,
          0.5096695247430424223356548135815582e-3,
          0.2231547584535793797614188036013401e-3,
          0.9945751278180853371459589003190170e-4,
          0.4492623673813314170020750240635786e-4,
          0.2050721277567069155316650397830591e-4,
          0.9439488275268395903987425104415055e-5,
          0.4374866789907487804181793223952411e-5,
          0.2039215753801366236781900709670839e-5,
          0.9551412130407419832857179772951265e-6,
          0.4492469198764566043294290331193655e-6,
          0.2120718480555466586923135901077628e-6,
          0.1004322482396809960872083050053344e-6,
          0.4769810169363980565760193417246730e-7,
          0.2271109460894316491031998116062124e-7,
          0.1083865921489695409107491757968159e-7,
          0.5183475041970046655121248647057669e-8,
          0.2483674543802478317185008663991718e-8,
          0.1192140140586091207442548202774640e-8,
          0.5731367241678862013330194857961011e-9,
          0.2759522885124233145178149692816341e-9,
          0.1330476437424448948149715720858008e-9,
          0.6422964563838100022082448087644648e-10,
          0.3104424774732227276239215783404066e-10,
          0.1502138408075414217093301048780668e-10,
          0.7275974480239079662504549924814047e-11,
          0.3527742476575915083615072228655483e-11,
          0.1711991790559617908601084114443031e-11,
          0.8315385841420284819798357793954418e-12,
          0.4042200525289440065536008957032895e-12,
          0.1966475631096616490411045679010286e-12,
          0.9573630387838555763782200936508615e-13,
          0.4664076026428374224576492565974577e-13,
          0.2273736960065972320633279596737272e-13,
          0.1109139947083452201658320007192334e-13 ##/* = (zeta(40+1)-1)/(40+1) */
          )

    c <- 0.2273736845824652515226821577978691e-12##/* zeta(N+2)-1 */

    r <- a
    ## if(abs(a) >= 0.5)
    a.lrg <- (abs(a) >= 0.5)
    r[a.lrg] <- lgamma(a[a.lrg] + 1)
    ## else {
    if(any(a.sml <- !a.lrg)) {
        a. <- a[a.sml]
        ## Abramowitz & Stegun 6.1.33 : for |x| < 2,
        ## <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * sum_{n=0..Inf} c_n (-x)^n
        ## where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
        ##
        ## Continued fraction convergence acceleration is used to compute
        ## lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
        lgam <- c * logcf(-a. / 2, N + 2, 1, tol_logcf)
        for (i in N:1)
            lgam <- coeffs[i] - a. * lgam

        ## return
        r[a.sml] <- (a. * lgam - eulers_const) * a. -
            log1pmx(a., tol_logcf = f.tol * tol_logcf, ...)
    }
    r
}## lgamma1p


lgamma1p_series <- function(x, k) {
    stopifnot(k == as.integer(k), 1 <= k, k <= 11)
    ## From Maple : ~/maple/gamma-asympt.txt

    ##  lG := series(ln(GAMMA(x+1)), x, 11);
    ##                   1    2  2   1          3    1    4  4
    ##  lG := -gamma x + -- Pi  x  - - Zeta(3) x  + --- Pi  x
    ##                   12          3              360

    ##       1          5    1     6  6   1          7     1     8  8
    ##     - - Zeta(5) x  + ---- Pi  x  - - Zeta(7) x  + ----- Pi  x
    ##       5              5670          7              75600

    ##       1          9     1      10  10    / 11\
    ##     - - Zeta(9) x  + ------ Pi   x   + O\x  /
    ##       9              935550

    ## Maple> hG := convert(convert(lG, polynom), horner);
    ## gives
    ## (-gamma +
    ##   (Pi*Pi/12 +
    ##    (-Zeta(3)/3+
    ##     (Pi*Pi*Pi*Pi/360+
    ##      (-Zeta(5)/5+
    ##       (Pi*Pi*Pi*Pi*Pi*Pi/5670+
    ##        (-Zeta(7)/7+
    ##         (Pi*Pi*Pi*Pi*Pi*Pi*Pi*Pi/75600+
    ##          (-Zeta(9)/9 +
    ##           Pi^10/935550 * x) * x) * x) * x) * x) * x) * x) * x) * x) * x

    useM <- (inherits(x, "mpfr"))
    if(useM) {
        stopifnot(requireNamespace("Rmpfr"))
        prec <- max(Rmpfr::getPrec(x))
        mpfr <- Rmpfr::mpfr
        zeta <- Rmpfr::zeta
    }
    ## "Euler's constant" gamma
    gamma <- if(useM) Rmpfr::Const("gamma", prec) else 0.57721566490153286
    if(useM) pi <- Rmpfr::Const("pi", prec) # else use R's pi # 3.1415926535897932

    if(k >= 2) {
        px <- (pi*x)^2 # = pi^2 x^2
        px4 <- px/4    # = pi^2 x^2 / 4
        if(k >= 3) {
            z3 <- if(useM) zeta(mpfr(3, prec)) else 1.202056903159594285
            x2 <- x*x
            if(k >= 5) {
                z5 <- if(useM) zeta(mpfr(5, prec)) else 1.036927755143369927
                if(k >= 6) {
                    I <- if(useM) mpfr(1, prec) else 1 # use I/3 = 1/3 in high precision
                    if(k >= 7) {
                        z7 <- if(useM) zeta(mpfr(7, prec)) else 1.008349277381922827
                        ## if(k >= 9) {
                        ##     z9 <- if(useM) zeta(mpfr(9, prec)) else 1.002008392826082214
                        ##     if(k >= 11) {
                        ##         z11 <- if(useM) zeta(mpfr(11, prec)) else 1.000494188604119464
                        ##     }
                        ## }

                    }
                }
            }
        }
    }
   ## return
    switch(k,
           -gamma * x,
           -gamma * x + px4/3,    # k = 2
           -gamma * x + px4/3 - x*x2/3*z3,    # k = 3
           -gamma * x - x*x2*(z3/3) + px4/3*(1 + px4/7.5),    # k = 4
           ##   pi^4/360 = (pi^2/4)^2 / 3 / 7.5  (3 * 16 * 7.5 == 360)
           ## k = 5 :
           -gamma * x - x*x2*(z3/3 + x2*(z5/5)) + px4/3*(1 + px4/7.5),
           ## k = 6 : (fractions: multiple of 2^-k are exact also in double prec)
           ##         360/16 = 22.5   5670/16 = 354.375
           -gamma * x - x*x2*(z3/3 + x2*(z5/5))            + px4*(1/3 + px4*(I/22.5 + px/354.375)),
           -gamma * x - x*x2*(z3/3 + x2*(z5/5+ x2*(z7/7))) + px4*(1/3 + px4*(I/22.5 + px/354.375)), # k = 7
           stop("Currently, only k <= 7  is implemented, but k=",k))
}

algdiv <- function(a,b) .Call(C_R_algdiv, as.double(a), as.double(b))



qnormAppr <- function(p) {
##' original (getting inaccurate when p --> 1)
  ## qnorm: normal quantile approximation for qnorm(p),  for p > 1/2
  ## -- to be used in  qbeta(.)
  ## The relative error of this approximation is quite ASYMMETRIC: mainly < 0
  C1 <- 2.30753
  C2 <- 0.27061
  C3 <- 0.99229
  C4 <- 0.04481
  r <- sqrt(-2*log(1-p))
  r - (C1 + C2 * r) / (1 + (C3 + C4 * r) * r)
}

qnormUappr <- function(p,
                       lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
                                        # ~= log(1-p) -- independent of lower.tail, log.p
                       lower.tail=FALSE, log.p=FALSE)
{
    ## qnorm: normal quantile approximation;
    ## -- to be used in  qbeta(.)
    ## The relative error of this approximation is quite ASYMMETRIC: mainly < 0
    C1 <- 2.30753
    C2 <- 0.27061
    C3 <- 0.99229
    C4 <- 0.04481

    if(missing(p)) { ## log.p unused;  lp = log(1-p)  <==>  e^lp = 1-p  <==>  p = 1 - e^lp
        p. <- -expm1(lp)
        ## swap p <--> 1-p -- so we are where approximation is better
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector

    } else {
        p. <- .D_qIv(p, log.p)
        ## swap p <--> 1-p -- so we are where approximation is better
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector
        p[swap] <- if(log.p) log1mexp(-p[swap]) else 1 - p[swap]
    }
    R <- r <- sqrt(-2 * lp)
    r.ok <- r < 1e10 ## e.g., for p == 1, r = Inf would give R = NaN
    ## for r >= 1e10:  R = r numerically in formula below
    r <- r[r.ok]
    R[r.ok] <- r - (C1 + C2 * r) / (1 + (C3 + C4 * r) * r)
    R[swap] <- -R[swap]
    R[p. == 1/2] <- 0
    R
}

qbetaAppr.1 <- function(a, p, q, lower.tail=TRUE, log.p=FALSE,
                        y = qnormUappr(a, lower.tail=lower.tail, log.p=log.p))
{
  ## Purpose: Approximate  qbeta(a, p,q) -- Abramowitz & Stegun (26.5.22)
  ##          qbeta(.) takes this only when  p>1 & q>1 : "Carter(1947), see AS 109, remark '5.'"
  ## -------------------------------------------------------------------------
  ## Arguments: a: percentage point, only used via qnorm*(a,..);  (p,q): beta parameters

  r <- (y * y - 3) / 6 ## lambda
  s <- 1 / (p + p - 1)
  t <- 1 / (q + q - 1)
  h <- 2 / (s + t)
  w <- y * sqrt(h + r) / h - (t - s) * (r + 5 / 6 - 2 / (3 * h))
  p / (p + q * exp(w + w))
}

qbetaAppr.3 <- function(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
{
   ##  a=alpha
    ## Purpose: Approximate  qbeta(a, p,q) -- for small a [log.p=TRUE: a far negative]
    ##  Inversion of   I_x(a,b) ~= x^a / (a B(a,b))

    ## CARE:  log(p) + logbeta == log(p * beta(p,q)) suffers from cancellation
    ## ----   for small p
    ## ---> look at Qab() above or *rather* also see experiments c12pBeta() , p.err.pBeta()
    ## in  ../tests/qbeta-dist.R
    ##              ============

    ## log.a = log("a") = log(pr_{lower_tail})
    log.a <- if(lower.tail) .D_log(a, log.p=log.p) else .D_LExp(a, log.p=log.p)

    pmin(1, exp(log.a +  log(p) + logbeta) / p)
}

qbetaAppr.2 <- function(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
{
    ## Purpose: Approximate  qbeta(a, p,q)
    ## log(1-"a") = log(pr_{upper_tail}):
    l1ma <- if(lower.tail) .D_LExp(a, log.p=log.p) else .D_log(a, log.p=log.p)
    ## pmax(0, ....) -- needed when  (l1ma + log(q) + logbeta) / q > 0, e.g., for
    ## e.g. qbetaAppr.2(1/4, 1/2, 1/2) gives -0.388
    -expm1((l1ma + log(q) + logbeta) / q)
}

qbetaAppr.4 <- function(a, p, q, lower.tail=TRUE, log.p=FALSE,
                        y = qnormUappr(a, lower.tail=lower.tail, log.p=log.p),
                        verbose=getOption("verbose"))
{
    ## Purpose: Approximate  qbeta(a, p,q); 'a' is only used via qnorm(a,..)
    r <- q + q
    t <- 1 / (9 * q)
    t <- r * (1 - t + y * sqrt(t))^3 # = \chi^2_alpha .. which "must be > 0", but I find: no!
    if(verbose) cat("chi^2[alpha] =: t = ",format(t), if(t <= 0)"t < 0 !!","\n")
    t <- (4 * p + r - 2) / t
    if(verbose) cat("t = (1 + x0)/(1 - x0) = ",format(t), if(t <= 1) "t <= 1 !!","\n")
    1 - 2 / (t + 1)
}

qbetaAppr.5 <- function(a, p, q, lower.tail=TRUE, log.p=FALSE,
                        y = qnormUappr(a, lower.tail=lower.tail, log.p=log.p),
                        logbeta = lbeta(p,q),
                        verbose=getOption("verbose"))
{
    ## Purpose: Approximate  qbeta(a, p,q) for p << q , really from  p --> 0 expansion
    ##
    ## based on future good  laBeta() ==robust== log(p * B(p,q))
    ##
    ## see MM's notes (on paper!)
}


qbetaAppr <- function(a, p, q, lower.tail=TRUE, log.p=FALSE,
                      y = qnormUappr(a, lower.tail=lower.tail, log.p=log.p),
                      logbeta = lbeta(p,q),
                      verbose = getOption("verbose") && length(a) == 1)
{
    ## Purpose: Approximate  qbeta(a, p,q) --- for  a <= 1/2
    ## -------------------------------------------------------------------------
    ## Arguments:
    ## -------------------------------------------------------------------------
    ## Author: After qbeta.c in  R  (--> AS 109...)
    if(length(p)!=1 || length(q)!=1)
        stop("'p' & 'q' must have length 1 !!")
    if(any(a > 1/2)) warning("a[.] > 1/2 (not thought for!)")

    if (p > 1 && q > 1) { ## Abramowitz & Stegun(26.5.22)
        if(verbose) cat("p,q > 1: Using qbetaAppr.1(): Abramowitz-Stegun (26.5.22)\n")
        ## 'a' is not needed here, just 'y' is
        qbetaAppr.1(,p,q, lower.tail=lower.tail, log.p=log.p, y=y)
    } else {
        if(verbose) cat("p or q <= 1 : ")
        r <- q + q
        ## more "robust" version of   t <- 1 / (9 * q) ; t <- r * (1 - t + y * sqrt(t))^3
        st <- 1/(3 * sqrt(q))
        t <- r * (1 - st*(st + y))^3

        ## NOTE: length(t) == length(y) == length(a)  is fulfilled

        ## vectorized in t { ==> logical vectors instead of simple if(..) .. else if(..) .. else ..
        ans <- t
        if(any(neg <- t <= 0)) { ## forget y(a) and t(q, y) from above; qbetaAppr.2():
            if(verbose)   cat("t <= 0:  using qbetaAppr.2()\n")
            ## log(1-"a") = log(pr_{upper_tail}):
            l1ma <- if(lower.tail) .D_LExp(a[neg], log.p=log.p) else .D_log(a[neg], log.p=log.p)
            ## FIXME -- still ? -- : "catastrophic"  in -expm1(xx) when xx > 0 ==> -expm1(xx) = 1-e^xx < 0 "nonsense"
            ans[neg] <- - expm1((l1ma + log(q) + logbeta) / q)
        }
        if(any(!neg)) { ##else
            t <- (4 * p + r - 2) / t
            if (any(L1 <- !neg & t <= 1)) {
                if(verbose) cat("t > 0; t' <= 1: qbetaAppr.3()\n")
                ## ans[L1] <- exp((log(a[L1] * p) + logbeta) / p)
                ans[L1] <- qbetaAppr.3(a[L1], p, q, lower.tail=lower.tail, log.p=log.p, logbeta=logbeta)
                ## ==> linear relationship on log-log scale:
                ##   log(ans) = 1/p * log(a) + (log p + log beta(p,q))/p
            }
            if (any(L2 <- !neg & t > 1)) { # the other case
                if(verbose) cat("t > 0; t' > 1 : qbetaAppr.4() = 1-2/(t'+1)\n")
                ans[L2] <- 1 - 2 / (t[L2] + 1)
            }
        }
        ans
    }
}

## Hmm, this is (sometimes!) not quite equivalent to C's :
## ex.  qbeta.R(0.5078, .01, 5) -> 2.77558e-15 but qbeta() now correctly gives 4.651e-31
qbeta.R	 <-  function(alpha, p, q,
                      lower.tail = TRUE, log.p = FALSE,
		      logbeta = lbeta(p,q),
		      low.bnd = 3e-308, up.bnd = 1-2.22e-16,
                      method = c("AS109", "Newton-log"),
                      tol.outer = 1e-15,
                      ## FIXME: (a,p,q) : and then uses (a, pp) .. hmm
		      f.acu = function(a,p,q) max(1e-300, 10^(-13- 2.5/pp^2 - .5/a^2)),
		      fpu = .Machine$ double.xmin,
		      qnormU.fun = function(u, lu) qnormUappr(p=u, lp=lu)
		    , R.pre.2014 = FALSE
		    , verbose = getOption("verbose") ## FALSE, TRUE, or 0, 1, 2, ..
		    , non.finite.report = verbose
		      )
{
### ----- following the C code in  ~/R/D/r-devel/R/src/nmath/qbeta.c
###
### ----------> make use of 'log.p' in	qnormUappr(* , log.p)	as well!

    ## Purpose:	 qbeta(.) in R -- edited;  verbose output
    ## -------------------------------------------------------------------------
    ## Arguments: alpha: percentage;  (p,q) Beta parameter
    ##	  low.bnd, up.bnd: algorithm parameter
    ##	acu, fpu:	accuracy;     .Machine$ double.xmin = 2.22e-308 for IEEE
    ##		SAE below is the most negative decimal exponent which does not
    ##		cause an underflow; a value of -308 or thereabouts will often be
    ##		OK in double precision.
    ##		sae <- -37 ;  fpu <- 10 ^ sae
    ##	  qnormU.fun(p) ~= qnorm(1 - p) = qnorm(p, ..., lower.tail=FALSE)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Apr 1997, 12:55
    ## ------ edited C-code  qbeta.c  into  R-code

    if(length(alpha)>1) stop("only works for length(1) 'alpha'")

    if (low.bnd >= up.bnd || low.bnd < 0 || up.bnd > 1)
	stop("MUST  0 <= low.bnd < up.bnd <= 1")

    ## Test for admissibility of parameters
    if (p < 0 || q < 0 || alpha < 0. || alpha > 1.) stop("DOMAIN_ERROR")

    ##  R_Q_P01_boundaries(alpha, 0, 1) :=
    if(log.p) {
        if(alpha == -Inf) return(if(lower.tail) 1 else 0)
        if(alpha ==   0 ) return(if(lower.tail) 0 else 1)
    } else {
        if (alpha == 0 || alpha == 1) return(if(lower.tail) alpha else 1-alpha)
    }

    p. <- .DT_qIv(alpha, lower.tail, log.p=log.p) # = lower_tail prob (in any case)

    if(R.pre.2014) {
        if(log.p && (p. == 0. || p. == 1.))
	return(p.) ## better than NaN or infinite loop;

    } else { ## can and *must* do better:
        if(log.p && (p. == 0. || p. == 1.))
            warning("p. is 0 or 1")
    }

    ## inflection point of pbeta() (outside [0,1] when p < 1 < q or q < 1 < p)
    ## x.ip <- (p-1) / (p+q -2)
    ## y.ip <- pbeta(x.ip, p, q)

    ## change tail if necessary, such that  beta(a,p,q)	 with  0 = a <= 1/2
    ## la := log(a), but without numerical cancellation:
    if (alpha <= 1/2) {
	a <- p.
	la <- if(lower.tail) .D_log(alpha, log.p) else .D_LExp(alpha, log.p)
        pp <- p; qq <- q; swap.tail <- FALSE
    } else {
	a <- .DT_CIv(alpha, lower.tail, log.p)
	la <- if(lower.tail) .D_LExp(alpha, log.p) else .D_log(alpha, log.p)
        pp <- q; qq <- p; swap.tail <- TRUE
    }
    if(verbose) cat(if(swap.tail)"SWAP tail" else "no swap", "\n")
    stopifnot( all.equal(log(a), la) )

    acu <- f.acu(a,p,q)
    ## set the exponent of acu to -2r-2 for r digits of accuracy
    if(acu <= 0 || acu > .1) acu   <- 1.0e-32 # 32: 15 digits accu.

    if(verbose) cat("logbeta=", formatC(logbeta,digits=19),
		    "; acu=",formatC(acu, digits=6) ,"\n")

    ## calculate the initial approximation 'xinbta' [FIXME: never designed for log.p=TRUE]
    xinbta <- qbetaAppr(a, pp,qq, y = qnormU.fun(u=a, lu=la), logbeta=logbeta,
                         verbose = verbose)
    xinbta0 <- xinbta # in case we want get back to it
    if(verbose) cat("initial 'xinbta' =", formatC(xinbta0, digits=18))
    if(!is.finite(xinbta)) {
	warning("** qbeta.R(): qbetaAppr() was not finite!   xinbta := 0.5")
        xinbta <- 0.5
    }

    ## Solve for x by a modified Newton-Raphson method,
    ## using the function  pbeta(.)

    r <- 1 - pp
    t <- 1 - qq
    yprev <- 0
    adj <- 1
    ##not used: prev <- 1
    if (xinbta < low.bnd) { xinbta <- 0.5; if(verbose)cat(" -- low.bnd  -> xi= 0.5\n") } else
    if (xinbta >  up.bnd) { xinbta <- 0.5; if(verbose)cat(" -- up.bnd   -> xi= 0.5\n") } else
    if(verbose) cat("\n")

    ##-- for single precision -- from CORRECTED	 AS 109:
    ##R	 iex <- as.integer( max(-5/pp^2 - 1/a^2 - 13, sae))
    ##NN acu <- max(fpu, 10 ^ (-5/pp^2 - 1/a^2 - 13))

    finish <- FALSE
    o.it <- 0
    method <- match.arg(method)

    switch(method,
           "AS109" =
    while(!finish) { ##-- Outer Loop
	y <- pbeta(xinbta, pp, qq)
	if(verbose) cat("y=pbeta(x..)=", formatC(y, digits=15, width=18, flag='-'))
	if(!is.finite(y)) {
            if(non.finite.report)
	    cat("  y=pbeta(): not finite: 'DOMAIN_ERROR'\n",
		"  xinbta =", format01prec(xinbta, digits=18),"\n")
	    return(xinbta) ##browser()
	}
	y <- (y - a) * exp(logbeta + r * log(xinbta) + t * log1p(-xinbta))
	if(verbose) cat(" y.n=(y-a)*e^..=", formatC(y))
	if(!is.finite(y)) {
            if(non.finite.report)
	    cat("  y = (y-a)*exp(...) not finite:\n",
		"  xinbta =", format01prec(xinbta, digits=18),"\n")
 	    return(xinbta) ##browser()
	}

	if (y * yprev <= 0) prev <- max(fpu, abs(adj))
	g <- 1
	in.it <- 0
	if(verbose>=2) cat("\n	Inner loop: ")
	for(in.it in 1:1000) {
	    adj <- g * y
	    if (abs(adj) < prev) { #-- current adjustment < last one which gave sign change
		if(verbose>=2)cat("<")
		tx <- xinbta - adj
		if (0 <= tx && tx <= 1) {
		    if (prev <= acu || abs(y) <= acu) { ## goto L_converged
			finish <- TRUE; break
		    }
		    if (tx != 0 && tx != 1) break
		}
	    } else if(verbose>=2) cat(".")
	    if(!finish) g <- g / 3
	} ##-- end inner loop
	if(verbose>=2) cat(if(in.it>10) "\n" else " ")
	if(verbose) cat(" ", in.it,"it --> adj=g*y=", formatC(adj),"\n")
	if (abs(tx - xinbta) <= tol.outer * xinbta) break # goto L_converged;
	if(verbose >= 2)
            cat("--Outer Loop: |tx - xinbta| > eps; tx=",formatC(tx, digits=18),
                " tx-xinbta =", formatC(tx-xinbta),"\n")
	xinbta <- tx
	yprev <- y
	o.it <- o.it + 1
    }, ##- end outer loop {method "AS109"}

    "Newton-log" = #--------------------------------- new --------
        ## but really, I want  Newton on   log-*log* scale --
    while(!finish) { ##-- Outer Loop
	y <- pbeta(xinbta, pp, qq, log.p = TRUE)
	if(verbose) cat("y=pbeta(xi,*, log.p=TRUE)=",
                        formatC(y, digits=15, width=18, flag='-'))
	if(!is.finite(y)) {
	    cat("  y=pbeta(): not finite: 'DOMAIN_ERROR'\n",
		"  xinbta =",formatC(xinbta,digits=18),"\n")
	    return(xinbta) ##browser()
	}
        ## y := g(xinbta) / g'(xinbta)  where g(x) =  log(alpha) - log pbeta(x, ..)
        ##                              ==>  g'(x) = - dbeta(x,.) / pbeta(x,.)
	y <- (y - la) * exp(y + (logbeta + r * log(xinbta) + t * log1p(-xinbta)))
	if(verbose) cat(" y.n=(y-a)*e^..=", formatC(y))
	if(!is.finite(y)) {
	    cat("  y = (y-a)*exp(...) not finite:\n",
		"  xinbta =",formatC(xinbta,digits=18),"\n")
 	    return(xinbta) ##browser()
	}

## 	if (y * yprev <= 0) prev <- max(fpu, abs(adj))
## 	g <- 1
## 	in.it <- 0
## 	if(verbose>=2) cat("\n	Inner loop: ")
## 	for(in.it in 1:1000) {
## 	    adj <- g * y
## 	    if (abs(adj) < prev) { #-- current adjustment < last one which gave sign change
## 		if(verbose>=2)cat("<")
## 		tx <- xinbta - adj
## 		if (0 <= tx && tx <= 1) {
## 		    if (prev <= acu || abs(y) <= acu) { ## goto L_converged
## 			finish <- TRUE; break
## 		    }
## 		    if (tx != 0 && tx != 1) break
## 		}
## 	    } else if(verbose>=2) cat(".")
## 	    if(!finish) g <- g / 3
## 	} ##-- end inner loop
## 	if(verbose>=2) cat(if(in.it>10)"\n" else " ")
## 	if(verbose) cat(" ",in.it,"it --> adj=g*y=",formatC(adj),"\n")
 	if(verbose) cat("\n")
        tx <- xinbta - y
	if (abs(tx - xinbta) < tol.outer * xinbta) break # goto L_converged;
	if(verbose>=2) cat("--Outer Loop: |tx - xinbta| > eps; tx=",formatC(tx, digits=18),
                           " tx-xinbta =", formatC(tx-xinbta),"\n")
	xinbta <- tx
## 	yprev <- y
	o.it <- o.it + 1
    }, ##- end outer loop, method = "Newton-log"
    ## otherwise:
    stop("unknown method ", dQuote(method)))

    ## L_converged:
    if(verbose) cat(" ", o.it,"outer iterations\n")
    if (swap.tail) 1 - xinbta else xinbta
}

