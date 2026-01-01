#### Bessel Function J_nu(x)
#### ----------------------------------------------


## It is possible to define the function by its Taylor series expansion around x = 0:

##  J_\alpha(x) = \sum_{m=0}^\infty \frac{(-1)^m}{m! \, \Gamma(m+\alpha+1)} {\left(\frac{x}{2}\right)}^{2m+\alpha}

## besselJ() - Definition as infinite sum  -- working for "mpfr" numbers, too!
## ---------
besselJs <-
    function(x, nu, nterm = 800, log = FALSE,
	     Ceps = if(isNum) 8e-16 else 2^(- x@.Data[[1]]@prec))
{
    ## Purpose: besselJ() primitively
    ## ----------------------------------------------------------------------
    ## Arguments: (x,nu) as besselJ;  nterm: number of terms for "infinite" sum
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: Dec 2014
    if(length(nu) > 1)
        stop(" 'nu' must be scalar (length 1)!")
    n <- length(x)
    if(n == 0) return(x)
    isNum <- is.numeric(x) || is.complex(x)
    if (nu < 0) {
	## Using Abramowitz & Stegun  9.1.2
	## this may not be quite optimal (CPU and accuracy wise)
        na <- floor(nu)
        return(if(!log)
                   (if(nu - na == 0.5) 0 else besselJs(x, -nu, nterm=nterm, Ceps=Ceps) * cospi(nu)) +
                   (if(nu      == na ) 0 else besselY (x, -nu                        ) * sinpi(nu))
               ## TODO: besselYs() series
               ## (if(nu      == na ) 0 else besselYs(x, -nu, nterm=nterm, Ceps=Ceps) * sinpi(nu))
               else ## same on log scale  ==> need lsum() ?
                   stop("besselJs(*, nu < 0, log = TRUE)  not yet implemented")
               )
    }

    j <- (nterm-1):0 # sum smallest first!
    sgns <- rep_len(if(nterm %% 2 == 0) c(-1,1) else c(1,-1), nterm)
    has0 <- any(i0 <- x == 0)
    x. <- if(has0) x[!i0] else x
    if(is.numeric(x) && is(nu, "mpfr")) {
	x. <- mpfr(x., precBits = max(64, .getPrec(nu)))
	isNum <- FALSE
    }
    l.s.j <- outer(j, x./2, function(X,Y) X*2*log(Y))##-> {nterm x n} matrix
    ##
    ## improve accuracy for lgamma(j+1)  for "mpfr" numbers
    ## -- this is very important [evidence: e.g. besselJs(10000, 1)]
    if(is(l.s.j, "mpfr"))
	j <- mpfr(j, precBits = max(sapply(l.s.j@.Data, slot, "prec")))
    else if(!isNum) j <- as(j, class(x))

    ## underflow (64-bit AMD) for x > 745.1332
    ## for large x, this already overflows to Inf :
    ## s.j <- outer(j, (x^2/4), function(X,Y) Y^X) ##-> {nterm x n} matrix
    ## s.j <- s.j / (gamma(j+1) * gamma(nu+1 + j))  but without overflow :
    log.s.j <- l.s.j - lgamma(j+1) - lgamma(nu+1 + j)
    s.j <-
        if(log) # NB: lsum() works on whole matrix
            lssum(log.s.j, signs=sgns) # == log(sum_{j} exp(log.s.j) )
        else ## log J_nu(x) -- trying to avoid overflow/underflow for large x OR large nu
            ## log(s.j) ; e..x <- exp(-x) # subnormal for x > 1024*log(2) ~ 710;
            exp(log.s.j)
    if(log) {
	if(any(lrgS <- log.s.j[1,] > log(Ceps) + s.j))
	    lapply(x.[lrgS], function(x)
		warning(gettextf(" 'nterm=%d' may be too small for x=%g", nterm, x),
			domain=NA))
	if(has0) {
	    sj <- x
	    sj[!i0] <- s.j
	    s.j <- sj
	}
	nu*log(x/2) + s.j
    } else { ## !log
	s <- colSums(sgns * s.j)
	if(!all(iFin <- is.finite(s)))
	    stop(sprintf("infinite s for x=%g", x.[!iFin][1]))
	if(any(lrgS <- s.j[1,] > Ceps * s))
	    lapply(x.[lrgS], function(x)
		warning(gettextf(" 'nterm=%d' may be too small for x=%g", nterm, x),
			domain=NA))
	if(has0) {
	    sj <- x
	    sj[!i0] <- s
	    s <- sj
	}
	(x/2)^nu * s
    }
}

##------------- Small nu  case  -- R-devel 2025-12 -----------------
## MM: see ~/R/MM/NUMERICS/Bessel/besselJ_small_nu.md
##         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bessJsmlNu.1 <- function(x, nu, nTrms, m.max) {
    stopifnot(length(x) == 1, length(nu) == 1, ## --> to be vectorized to bJsmlNu()
              nTrms == as.integer(nTrms), length(nTrms) == 1, 0 <= nTrms, nTrms <= 2,
              m.max == as.integer(m.max), length(m.max) == 1, 1 <= m.max)
    J0 <- besselJ(x, 0)
    if(nTrms == 0)
        return(J0)
    ## else nTrms >= 1 ,  i.e., (for now) in { 1 , 2 }:
    m <- 0:m.max; sig.m <- rep_len(c(1, -1), m.max+1) # = (-1)^m
    x2 <- x/2 # ldexp(x, -1)
    m.fac <- factorial(m)
    a_m <- sig.m * x2^(2*m)/m.fac^2 # := (-1)^m (x/2)^{2m}/(m!)^2
    ## Compute the nu = nu^1 "coefficient" (an infinite series sum)
    ## s1 <- 0
    ##- for(m in 0:m.max) -->
    ##    s1 <- s1 + sum(a_m *  ...........)

    ## J1 :=  ν · Σ_{m≥0} a_m ( ln(x/2) − ψ(m+1) )
    s1 <- sum( a_m * (log(x2) - digamma(m+1)))
    J1 <- nu * s1
    if(nTrms == 1)
        return(J0 + J1)
    ## else nTrms >= 2
    ## Compute the nu^2 "coefficient" (an infinite series sum)
    ## s2 <- 0
    ## for(m in 0:m.max)
    ##    s2 <- s2 +  ........

    ## J2 :=  ν^2 · Σ_{m≥0} a_m · ( 1/2 ln^2(x/2) − ln(x/2) ψ(m+1) + 1/2(ψ(m+1)^2 − ψ'(m+1)) )
    s2 <- sum( a_m * (log(x2)^2/2 - log(x2)* digamma(m+1) + 1/2*(digamma(m+1)^2 - trigamma(m+1))) )
    J2 <- nu^2 * s2
    if(nTrms > 2)
        stop(" 'nTrms' must be in {0, 1, 2} currently")
    ## return
    J0 + J1 + J2
}

bessJsmlNu <- Vectorize(bessJsmlNu.1, c("x", "nu"))
