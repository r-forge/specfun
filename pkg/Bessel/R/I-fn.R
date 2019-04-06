###--- This used to be part /u/maechler/R/MM/NUMERICS/bessel-fn.R

#### "Modified" Bessel Function I_nu(x)
####  ----------------------------------------------
## besselI() - Definition as infinite sum  -- working for "mpfr" numbers, too!
## ---------
besselIs <-
    function(x, nu, nterm = 800, expon.scaled = FALSE, log = FALSE,
	     Ceps = if(isNum) 8e-16 else 2^(- x@.Data[[1]]@prec))
{
    ## Purpose: besselI() primitively
    ## ----------------------------------------------------------------------
    ## Arguments: (x,nu) as besselI;  nterm: number of terms for "infinite" sum
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 21 Oct 2002, 21:59
    if(length(nu) > 1)
        stop(" 'nu' must be scalar (length 1)!")
    n <- length(x)
    if(n == 0) return(x)
    j <- (nterm-1):0 # sum smallest first!
    if(is.numeric(x) && is(nu, "mpfr")) {
	x <- mpfr(x, precBits = max(64, .getPrec(nu)))
	isNum <- FALSE
    } else
	isNum <- is.numeric(x) || is.complex(x)
    l.s.j <- outer(j, (x/2), function(X,Y) X*2*log(Y))##-> {nterm x n} matrix
    ##
    ## improve accuracy for lgamma(j+1)  for "mpfr" numbers
    ## -- this is very important [evidence: e.g. bI(10000, 1)]
    if(is(l.s.j, "mpfr"))
	j <- mpfr(j, precBits = max(sapply(l.s.j@.Data, slot, "prec")))
    else if(!isNum) j <- as(j, class(x))

    ## underflow (64-bit AMD) for x > 745.1332
    ## for large x, this already overflows to Inf :
    ## s.j <- outer(j, (x^2/4), function(X,Y) Y^X) ##-> {nterm x n} matrix
    ## s.j <- s.j / (gamma(j+1) * gamma(nu+1 + j))  but without overflow :
    log.s.j <-
	if(expon.scaled)
	    l.s.j - rep(abs(Re(x)), each=nterm) - lgamma(j+1) - lgamma(nu+1 + j)
	else
	    l.s.j	 	                - lgamma(j+1) - lgamma(nu+1 + j)
    if(log) {
	s.j <- lsum(log.s.j) # == log(sum_{j} exp(log.s.j) ); NB: lsum() works on whole matrix
	if(any(i0 <- x == 0)) s.j[i0] <- 0
	if(any(lrgS <- Re(log.s.j[1,]) > Re(log(Ceps) + s.j)))
	    lapply(x[lrgS], function(x)
		warning(gettextf(" 'nterm=%d' may be too small for %s", nterm,
				 paste(format(x), collapse=", ")),
			domain=NA))
	nu*log(x/2) + s.j
    } else { ## !log
	s.j <- ## log I_nu(x) -- trying to avoid overflow/underflow for large x OR large nu
	    ## log(s.j) ; e..x <- exp(-x) # subnormal for x > 1024*log(2) ~ 710;
	    exp(log.s.j)
	s <- colSums(s.j)
	if(any(i0 <- x == 0)) s[i0] <- 0
	if(!all(iFin <- is.finite(s)))
	    stop(sprintf("infinite s for x=%g", x[!iFin][1]))
	if(any(lrgS <- Re(s.j[1,]) > Re(Ceps * s)))
	    lapply(x[lrgS], function(x)
		warning(gettextf(" 'nterm=%d' may be too small for %s", nterm,
				 paste(format(x), collapse=", ")),
			domain=NA))
	(x/2)^nu * s
    }
}

## old name [back compatibility]:
bI <- function(x, nu, nterm = 800, expon.scaled = FALSE, log = FALSE,
               Ceps = if(isNum) 8e-16 else 2^(- x@.Data[[1]]@prec))
{
    .Deprecated("besselIs")
    isNum <- is.numeric(x) || is.complex(x)
    besselIs(x, nu, nterm=nterm, expon.scaled=expon.scaled, log=log, Ceps=Ceps)
}

###--------------- besselI() for large x or also large nu ---
###
###  see --> ~/R/MM/NUMERICS/bessel-large-x.R  for code usage
###          ================================

besselIasym <- function(x, nu, k.max = 10, expon.scaled=FALSE, log=FALSE)
{
    ## Purpose: Asymptotic expansion of Bessel I_nu(x) function   x -> oo
    ##        by Abramowitz & Stegun (9.7.1), p.377 :
    ##
    ## I_a(z) = exp(z) / sqrt(2*pi*z) * f(z,..)  where
    ##   f(z,..) = 1 - (mu-1)/ (8*z) + (mu-1)(mu-9)/(2! (8z)^2) - ...
    ##           = 1- (mu-1)/(8z)*(1- (mu-9)/(2(8z))*(1-(mu-25)/(3(8z))*..))
    ## where  mu = 4*a^2  *and*  |arg(z)| < pi/2
    ## ----------------------------------------------------------------------
    ## Arguments: x, nu, expon.scaled:  as besselI()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr, 21 Nov 2008

    ## Note that for each x the series eventually *DIVERGES*
    ## it should be "stopped" as soon as it *converged* (:-)

    stopifnot(k.max == round(k.max))

    ## I_\nu(z) = exp(z) / sqrt(2*pi*z) * f(z, nu)

    ## First compute  f(x, nu) to order k.max
    ## f <- 1
    d <- 0 ## d = 1 - f  <==>  f = 1 - d
    if(k.max >= 1) {
	## m. <- 4*nu^2
	x8 <- 8*x
        for(k in k.max:1) {
            ## mu <- 4*nu^2; d <- (1 - d)*(mu - (2*k-1)^2)/(k*x8)
            ## mu - (2k-1)^2 = (2nu - (2k-1)) (2nu + (2k-1)) =
            ##               = (2(nu-k)+1)(2(nu+k)-1)
            d <- (1 - d)*((2*(nu-k)+1)*(2*(nu+k)-1))/(k*x8)
        }
    }
    if(expon.scaled)
        sx <- x -  abs(if(is.complex(x)) Re(x) else x)
    pi2 <- 2* (if(inherits(x, "mpfr")) Rmpfr::Const("pi", max(.getPrec(x)))
               else pi)
    if(log) {
        ## f = 1 - d  ==> log(f) = log1p(-d) :
           (if(expon.scaled) sx else x) + log1p(-d) - log(pi2*x) / 2
    } else {
        exp(if(expon.scaled) sx else x) *   (1-d)  / sqrt(pi2*x)
    }
}

besselKasym <- function(x, nu, k.max = 10, expon.scaled=FALSE, log=FALSE)
{
    ## Purpose: Asymptotic expansion of Bessel K_nu(x) function   x -> oo
    ##        by Abramowitz & Stegun (9.7.2), p.378 :
    ##
    ## K_nu(z) = exp(-z) * sqrt(pi/(2*z)) * f(z,..)  where
    ##   f(z,..) = 1 + (mu-1)/ (8*z) + (mu-1)(mu-9)/(2! (8z)^2) + ...
    ##           = 1 + (mu-1)/(8z)*(1 + (mu-9)/(2(8z))*(1 + (mu-25)/(3(8z))*..))
    ## where  mu = 4*nu^2  *and*  |arg(z)| < pi/2
    ## ----------------------------------------------------------------------
    ## Arguments: x, nu, expon.scaled:  as besselK()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 15. Dec 2018

    ##__NB: This *very* similar to besselIasym():  1) different initial factor; 2) *not* alternating

    ## Note that for each x the series eventually *DIVERGES*
    ## it should be "stopped" as soon as it *converged* (:-)

    stopifnot(k.max == round(k.max))

    ## K_\nu(z) = exp(-z) * sqrt(pi/(2*z)) * f(z, nu)

    ## First compute  f(x, nu) to order k.max
    ## f <- 1 -- via  d := f - 1  <==>  f = 1 + d
    d <- 0
    if(k.max >= 1) {
	x8 <- 8*x
        for(k in k.max:1) {
            ## mu <- 4*nu^2; d <- (1 + d)*(mu - (2*k-1)^2)/(k*x8)
            ## mu - (2k-1)^2 = (2nu - (2k-1)) (2nu + (2k-1)) =
            ##               = (2(nu-k)+1)(2(nu+k)-1)
            d <- (1 + d)*((2*(nu-k)+1)*(2*(nu+k)-1))/(k*x8)
        }
    }
    pi.2 <- (if(inherits(x, "mpfr")) Rmpfr::Const("pi", max(.getPrec(x))) else pi)/2
    if(log) {
        ## f = 1 + d  ==> log(f) = log1p(d) :
        (if(expon.scaled) 0 else   -x   ) + log1p(d)  + (log(pi.2) - log(x)) / 2
    } else {
        (if(expon.scaled) 1 else exp(-x)) *  (1+d)    * sqrt(pi.2/x)
    }
}

besselI.ftrms <- function(x, nu, K = 20)
{
    ## Purpose: all *Terms* in besselIasym()
    ##        by Abramowitz & Stegun (9.7.1), p.377 :
    ##   f(x,..) = 1 - (mu-1)/ (8*x) + (mu-1)(mu-9)/(2! (8x)^2) - ...
    ##           = 1- (mu-1)/(8x)*(1- (mu-9)/(2(8x))*(1-(mu-25)/(3(8x))*..))
    ## where  mu = 4*nu^2
    ## ----------------------------------------------------------------------
    ## Arguments: x, nu:  as besselI()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr, 21 Nov 2008

    stopifnot(length(x) == 1, x > 0, length(nu) == 1,
              length(K) == 1, K == round(K), K >= 0)

    ## I_\nu(x) = exp(x) / sqrt(2*pi*x) * f(x, nu)

    kk <- seq_len(K)
    mu <- 4*nu^2
    x8 <- 8*x
    multf <- - (mu - (2*kk-1)^2)/(kk*x8)

    ## ser <- cumprod(multf)
    ## now, for k-term approximation f_k  of  f(x,nu)   we have
    ##  f_k = 1 + sum(ser[1:k]), i.e.    ##  (f_k)_k = 1 + cumsum(c(0, ser))[k+1]  for k= 0,1,...,K

    ## ser :=
    cumprod(multf)
}


###----------------------------------------------------------------------------

### When  BOTH  nu and x  are large :
##        ^^^^ (in practice, works fine already in some cases of |x| large)

besselI.nuAsym <- function(x, nu, k.max, expon.scaled=FALSE, log=FALSE)
{
    ## Purpose: Asymptotic expansion of Bessel I_nu(x) function
    ##        when BOTH  nu and x  are large

    ## Abramowitz & Stegun , p.378, __ 9.7.7. __
    ##
    ## I_nu(nu * z) ~ 1/sqrt(2*pi*nu) * exp(nu*eta)/(1+z^2)^(1/4) *
    ##                * {1 + u_1(t)/nu + u_2(t)/nu^2 + ... }
    ## where
    ## __ 9.7.11 __
    ##	       t := 1 / sqrt(1 + z^2)  = 1/sz
    ##       eta := sqrt(1 + z^2) + log(z / (1 + sqrt(1+z^2))) = sz + log(z / (1 + sz))
    ## with   sz := sqrt(1 + z^2)
    ##
    ## and u_k(t)  from  p.366  __ 9.3.9 __
    ##
    ## u0(t) = 1
    ## u1(t) = (3*t - 5*t^3)/24
    ## u2(t) = (81*t^2 - 462*t^4 + 385*t^6)/1152
    ## ... up to  u4(t)

    ## with recursion  9.3.10    for  k = 0, 1, .... :
    ##
    ## u_{k+1}(t) = t^2/2 * (1 - t^2) * u'_k(t) +
    ##            1/8  \int_0^t (1 - 5*s^2)* u_k(s) ds

    ## ----------------------------------------------------------------------
    ## Arguments: x, nu, expon.scaled:  as besselI()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 22 Nov 2008, 15:22

    stopifnot(k.max == round(k.max),
              0 <= k.max, k.max <= 5)

    z <- x/nu # -> computing I_nu(z * nu)
    sz <- sqrt(1 + z^2) ## << "FIXME": use hypot(.)
    t <- 1/sz
    ## if(expon.scaled) ## scale by   * exp(-abs(Re(x))):
    ##     sx <- x -  abs(if(is.complex(x)) Re(x) else x)
    eta <- (if(expon.scaled) { ## <==> * exp(-|Re(x)|) ==  exp(- |Re(z) * nu|)  <==> "eta - |Re(z)|"
		## For real z, have
		##  sz - |z| = sqrt(1 + z^2) - |z|  =!=  1/(sqrt(..) + |z|) = 1/(sz + |z|);
		## in complex case sz - |x| =  w/(sz + |x|)  with w := 1-y^2 + 2xy i; z = x + iy
		w <- if(is.complex(z)) {
			 x. <- Re(z); y <- Im(z)
			 1-y^2 + 2i*x.*y # complex(real = 1-y^2, imaginary = 2*x.*y)
		     } else 1
		w/(sz + abs(if(is.complex(z)) x. else z))
	    }
            else sz) + log(z / (1 + sz))
    ## I_nu(nu * z) ~ 1/sqrt(2*pi*nu) * exp(nu*eta)/(1+z^2)^(1/4) *
    ##                * {1 + u_1(t)/nu + u_2(t)/nu^2 + ... }
    if(k.max == 0)
        d <- 0
    else { ## k.max >= 1 --- Find the  Debye  polynomials u_j(t) in ../misc/QRMlib_src_bessel.c
           ## or in a newer GSL, but much more hidden ..../gsl-2.5/specfunc/debye.c
        t2 <- t^2
        u1.t <- (t*(3 - 5*t2))/24
        d <-
            if(k.max == 1) {
                u1.t/nu
            } else { ## k.max >= 2
                u2.t <- ##(81*t^2 - 462*t^4 + 385*t^6)/1152
                    t2*(81 + t2*(-462 + t2*385)) / 1152
                if(k.max == 2)
                    (u1.t + u2.t/nu)/nu
                else { ## k.max >= 3
                    u3.t <- t*t2*(30375 +
                                  t2*(-369603 +
                                      t2*(765765 - t2*425425)))/ 414720
                    if(k.max == 3)
                        (u1.t + (u2.t + u3.t/nu)/nu)/nu
                    else { ## k.max >= 4
                        t4 <- t2*t2
			u4.t <- t4*(4465125 +
                                    t2*(-94121676 +
                                        t2*(349922430 +
                                            t2*(-446185740 + t2*185910725))))/39813120
			if(k.max == 4)
			    (u1.t + (u2.t + (u3.t + u4.t/nu)/nu)/nu)/nu
                        else { ## k.max >= 5
                            u5.t <- t*t4*(1519035525 +
                                          t2*(-49286948607 +
                                              t2*(284499769554 +
                                                  t2*(-614135872350 +
                                                      t2*(566098157625 - t2*188699385875))
                                              )))/6688604160
                            if(k.max == 5)
                                (u1.t + (u2.t + (u3.t + (u4.t + u5.t/nu)/nu)/nu)/nu)/nu
                            else
                                stop("k.max > 5: not yet implemented (but should NOT happen)")
                        }
                    }
                }
            }
    }

    pi2 <- 2* (if(inherits(x, "mpfr")) Rmpfr::Const("pi", max(.getPrec(x)))
               else pi)
    if(log) {
        log1p(d) + nu*eta - (log(sz) + log(pi2*nu))/2
    } else {
        (1+d) * exp(nu*eta) / sqrt(pi2*nu*sz)
    }
} ## {besselI.nuAsym}

besselK.nuAsym <- function(x, nu, k.max, expon.scaled=FALSE, log=FALSE)
{
    ## Purpose: Asymptotic expansion of Bessel K_nu(x) function
    ##        when BOTH  nu and x  are large

    ## Abramowitz & Stegun , p.378, __ 9.7.8. __
    ##
    ## K_nu(nu * z) ~ sqrt(pi/(2*nu)) * exp(-nu*eta)/(1+z^2)^(1/4) *
    ##                * {1 - u_1(t)/nu + u_2(t)/nu^2 - ... }

    ## { see besselI.nuAsym() above, for t, eta, u_k ...}
    ## ----------------------------------------------------------------------
    ## Arguments: x, nu, expon.scaled:  as besselK()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 23 Nov 2009

    stopifnot(k.max == round(k.max),
              0 <= k.max, k.max <= 5)

    z <- x/nu # -> computing K_nu(z * nu)
    sz <- sqrt(1 + z^2) ## << "FIXME": use hypot(.)
    t <- 1/sz
    eta <- (if(expon.scaled) ## <==> * exp(-x)  ==  exp(- z * nu)  <==> "eta - z"
            ## sz - z = sqrt(1 + z^2) - z  =  1/(sqrt(..) + z)
            1/(z + sz) else sz) +
                log(z / (1 + sz))

    ## K_nu(nu * z) ~  [see above]
    if(k.max == 0)
        d <- 0
    else { ## k.max >= 1 --- Find the  Debye  polynomials u_j(t) ..... [see above!]
        t2 <- t^2
        u1.t <- (t*(3 - 5*t2))/24
        ## NB: Difference here for K() to I() above is *alternating signs*: "-" for odd uj()
        d <-
            if(k.max == 1) {
                - u1.t/nu
            } else { ## k.max >= 2
                u2.t <- ##(81*t^2 - 462*t^4 + 385*t^6)/1152
                    t2*(81 + t2*(-462 + t2*385)) / 1152
                if(k.max == 2)
                    (- u1.t + u2.t/nu)/nu
                else { ## k.max >= 3
                    u3.t <- t*t2*(30375 +
                                  t2*(-369603 +
                                      t2*(765765 - t2*425425)))/ 414720
                    if(k.max == 3)
                        (- u1.t + (u2.t - u3.t/nu)/nu)/nu
                    else { ## k.max >= 4
                        t4 <- t2*t2
                        u4.t <- t4*(4465125 +
                                    t2*(-94121676 +
                                        t2*(349922430 +
                                            t2*(-446185740 + t2*185910725))))/39813120
                        if(k.max == 4)
                            (- u1.t + (u2.t + (-u3.t + u4.t/nu)/nu)/nu)/nu
                        else { ## k.max >= 5
                            u5.t <- t*t4*(1519035525 +
                                          t2*(-49286948607 +
                                              t2*(284499769554 +
                                                  t2*(-614135872350 +
                                                      t2*(566098157625 - t2*188699385875))
                                              )))/6688604160
                            if(k.max == 5)
                                (- u1.t + (u2.t + (-u3.t + (u4.t - u5.t/nu)/nu)/nu)/nu)/nu
                            else
                                stop("k.max > 5: not yet implemented (but should NOT happen)")
                        }
                    }
                }
            }
    }

    if(log) {
        log1p(d) - nu*eta - (log(sz) - log(pi/(2*nu)))/2
    } else {
        (1+d) * exp(-nu*eta)*sqrt(pi/(2*nu * sz))
    }
} ## {besselK.nuAsym}
