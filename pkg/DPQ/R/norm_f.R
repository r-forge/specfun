#### qnorm(), pnorm() etc
#### --------------------

### Bounds on  1-Phi(.) = pnorm(*, lower.tail=FALSE) -- typically generalized Mill's ratios

### From Lutz Dümbgen (2010)'s arXiv
invisible("https://arxiv.org/abs/1012.2063")

## Sampford(1953)'s upper bound for 1-Phi(x)  --- Dumbgen's (5), p.2 :
pnormU_S53 <- function(x, lower.tail=FALSE, log.p=FALSE) {
    stopifnot(x >= 0)
    ## upper tail and *not* log.p :
    ##  4*dnorm(x) / (sqrt(8+x^2) + 3*x)

    ## log.p=TRUE and upper tail, i.e.  !lower.tail :
    ## log(4) + dnorm(x, log=TRUE) - log(sqrt(8+x^2) + 3*x)
    ## with less quick overflow for very large 'x' :
    ## sqrt(8+x^2) + 3x = |x| sqrt(1 + 8/x^2) + 3x = x(sqrt(1 + 8/x^2) + 3)
    ##
    ## dnorm(x, log=TRUE) - log(x) + (log(4) - log(sqrt(1 + 8/x^2) + 3))
    r <- dnorm(x, log=TRUE) - log(x) + log(4 / (3 + sqrt(1 + (8/x)/x)))
    if(log.p) {
        if(lower.tail) ## log(1 - exp(r)) = log1mexp(-r)
            log1mexp(-r)
        else ## upper tail: log(1 - (1 - exp(r))) = r
            r
    } else {
        if(lower.tail) -expm1(r) else exp(r)
    }
}

## Duembgen's lower bound (10), p.6
## { which is strictly better than Komatu(1955)'s lower bound (3) }
pnormL_LD10 <- function(x, lower.tail=FALSE, log.p=FALSE) {
    stopifnot(x > 0)
    ## non-log, upper tail :
    ## 1-Phi(x) >  ~=~ pi*dnorm(x) / ((pi-1)*x + sqrt(2*pi + x^2))
    ## log.p=TRUE and upper tail, i.e.  !lower.tail :
    r <- dnorm(x, log=TRUE) - log(x) + log(pi / (pi + sqrt(1 + (2*pi/x)/x) -1))
    if(log.p) {
        if(lower.tail) ## log(1 - exp(r)) = log1mexp(-r)
            log1mexp(-r)
        else
            r
    } else {
        if(lower.tail) -expm1(r) else exp(r)
    }
}



### R version of   ~/R/D/r-devel/R/src/nmath/qnorm.c

## Mathlib : A C Library of Special Functions
## Copyright (C) 2000--2020 The R Core Team
## Copyright (C) 1998       Ross Ihaka
## based on AS 111 (C) 1977 Royal Statistical Society
## and   on AS 241 (C) 1988 Royal Statistical Society

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/

## SYNOPSIS

##     double qnorm5(double p, double mu, double sigma,
##     	      int lower_tail, int log_p)
##           {qnorm (..) is synonymous and preferred inside R}

## DESCRIPTION

##     Compute the quantile function for the normal distribution.

##     For small to moderate probabilities, algorithm referenced
##     below is used to obtain an initial approximation which is
##     polished with a final Newton step.

##     For very large arguments, an algorithm of Wichura is used.

## REFERENCE

##     Beasley, J. D. and S. G. Springer (1977).
##     Algorithm AS 111: The percentage points of the normal distribution,
##     Applied Statistics, 26, 118-121.

##     Wichura, M.J. (1988).
##     Algorithm AS 241: The Percentage Points of the Normal Distribution.
##     Applied Statistics, 37, 477-484.


qnormR1 <- function(p, mu=0, sd=1, lower.tail=TRUE, log.p=FALSE,
                    trace = 0) {
    stopifnot(length(p) == 1, length(mu) == 1, length(sd) == 1,
              length(lower.tail) == 1, length(log.p) == 1,
              !is.na(lower.tail), !is.na(log.p))
    if(is.na(p) || is.na(mu) || is.na(sd)) return(p+mu+sd)
    ## R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
    if(p == .D_0(log.p)) return(if(lower.tail) -Inf else  Inf)
    if(p == .D_1(log.p)) return(if(lower.tail)  Inf else -Inf)
    if(p < .D_0(log.p) ||
       p > .D_1(log.p)) { warning("p out of range"); return(NaN) }

    if(sd < 0) { warning("sd < 0"); return(NaN) }
    if(sd == 0)	return(mu)

    p. <- .DT_qIv(p, lower.tail, log.p=log.p) # = lower_tail prob (in any case)
    q = p. - 0.5;

    if(trace)
        cat(sprintf("qnormR1(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
                    p,mu,sd, lower.tail, log.p, q))

## -- use AS 241 --- */
##  double ppnd16_(double *p, long *ifault)*/
##       ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

##       Produces the normal deviate Z corresponding to a given lower
##       tail area of P; Z is accurate to about 1 part in 10**16.

##       (original fortran code used PARAMETER(..) for the coefficients
##        and provided hash codes for checking them...)

    if (abs(q) <= .425) { ## 0.075 <= p <= 0.925 */
        r  <- .180625 - q * q;
	val <-
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608
            ) /
            (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.)
    }
    else { ## closer than 0.075 from {0,1} boundary */
           ##   r := log(p~);  p~ = min(p, 1-p) < 0.075 :
	r <-
            if(log.p && ((lower.tail && q <= 0) || (!lower.tail && q > 0))) {
                p
            } else {
                log( if(q > 0) .DT_CIv(p, lower.tail, log.p) ## 1-p
                     else p. ) # p. = R_DT_Iv(p) ^=  p
            }
	###  r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 ) :
        r  <- sqrt(-r)
        if(trace)
            cat(sprintf("\t close to 0 or 1: r = %7g\n", r))
        if(is.na(r))
            warning("r = sqrt( - log(min(p,1-p)) )  is NA -- \"we have lost\"")

        if (!is.na(r) && r <= 5.) { ## <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r <- r + -1.6;
            val <- (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                  1.42343711074968357734
            ) / (((((((r *
                       1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                      r + .0151986665636164571966) * r +
                     .14810397642748007459) * r + .68976733498510000455) *
                   r + 1.6763848301838038494) * r +
                  2.05319162663775882187) * r + 1.)
        }
        else { ## p is very close to  0 or 1:  r > 5 <==> min(p,1-p) < exp(-25) = 1.3888..e-11
            r <- r + -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772
            ) / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

	if(q < 0)
	    val  <- -val;
        ## return (q >= 0.)? r : -r ;*/
    }
    return(mu + sd * val)
} # qnormR1()

qnormR  <- Vectorize(qnormR1, c("p", "mu", "sd"))
##==

