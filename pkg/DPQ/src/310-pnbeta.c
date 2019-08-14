/* Produced by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 */
/* 310-CumNonCentralBeta.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip

-- be brave, try without   #include "f2c.h" ------*/

#include <Rmath.h>
#include <R.h>

double ncbeta_(double a, double b, double lambda, double x, double errmax,
	       int *ifault);

double betanc_(double x, double a, double b, double lambda,
	       double errmax, int itrmax, int *ifault);

/* call this from R : */
void sub_ncbeta(double *a, double *b, double *lambda, double *x, double *errmax,
		int *ifault, double *res, int *n_res)
{
    int i;
    for(i=0; i < *n_res; i++)
	res[i] = ncbeta_(a[i], b[i], lambda[i], x[i],
			 *errmax, ifault);
}

double ncbeta_(double a, double b, double lambda, double x, double errmax,
	       int *ifault)
{
    /* Initialized data */

    static int itrmax = 100; /* only for AS 226 currently */

    static double zero = 0.;
    static double half = 0.5;
    static double one =  1.;
    static double five = 5.;

    /* System generated locals */
    double ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(double), log(double), exp(double), lgamma(double);
/*
      double precision function betain(x, p, q, beta, ifault)

      double precision function gammad(x, p, ifault)
*/

    /* Local variables */
    static double c__;
    static int i__, j, m;
    static double q, r__, s, t, s0, s1, t0, t1, fx, gx;
    static int xj;
    static double ebd, sum, beta, temp, psum;
    static int iter1, iter2;
    static double errbd, ftemp;
    static int iterhi, iterlo;


/*       ALGORITHM AS 310 APPL. STATIST. (1997), VOL. 46, NO. 1

       Computes the cumulative distribution function of a
       non-central beta random variable



       Local variable XJ gives the number of iterations taken */




    ret_val = x;

/*       Check for admissibility of parameters */

    *ifault = 3;
    if (lambda <= zero || a <= zero || b <= zero) {
	return ret_val;
    }
    *ifault = 2;
    if (x < zero || x > one) {
	return ret_val;
    }
    *ifault = 1;
    if (x == zero || x == one) {
	return ret_val;
    }
    *ifault = 0;

    c__ = lambda * half;
    xj = zero;

    if (lambda < 54.f) {

/*       AS 226 as it stands is sufficient in this situation */

	ret_val = betanc_(x, a, b, lambda, errmax, itrmax, ifault);
	return ret_val;
    } else {
	m = (int) (c__ + half);
	t = five * sqrt((double) m);
	iterlo = m - t;
	iterhi = m + t;
	t = -c__ + m * log(c__) - lgamma(m + one);
	q = exp(t);
	r__ = q;
	psum = q;
	r__1 = a + m;
	r__2 = a + m + b;
	beta = lgamma(a + m) + lgamma(b) - lgamma(a + m + b);
	s1 = (a + m) * log(x) + b * log1p( - x) - log(a + m) - beta;
	gx = exp(s1);
	fx = gx;
	r__1 = a + m;
/*	temp = betain_(x, &r__1, b, &beta, ifault); */
	temp = pbeta(x, a + m, b, /*lower= */TRUE, /*log_p=*/FALSE);
	ftemp = temp;
	xj += one;
	sum = q - temp;
	iter1 = m;

/*      The first set of iterations starts from M and goes downwards */

L20:
	if (iter1 < iterlo) {
	    goto L30;
	}
	if (q < errmax) {
	    goto L30;
	}
	q -= iter1 / c__;
	xj += one;
	gx = (a + iter1) / (x * (a + b + iter1 - one)) * gx;
	iter1 -= one;
	temp += gx;
	psum += q;
	sum += q * temp;
	goto L20;
L30:
	t0 = lgamma(a + b) - lgamma(a + one) - lgamma(b);
	s0 = a * log(x) + b * log1p( - x);

	for (i__ = 1; i__ <= iter1; ++i__) {
	    j = i__ - one;
	    s += exp(t0 + s0 + j * log(x));
	    t1 = log(a + b + j) - log(a + one + j) + t0;
	    t0 = t1;
	}

/*       Compute the first part of error bound */

	errbd = pgamma(c__, (double)iter1, 1.,/*lower= */FALSE, /*log_p=*/FALSE)
	    * (temp + s);
	q = r__;
	temp = ftemp;
	gx = fx;
	iter2 = m;
L50:
	ebd = errbd + (one - psum) * temp;
	if (ebd < errmax || iter2 >= iterhi) {
	    goto L60;
	}
	iter2 += one;
	xj += one;
	q = q * c__ / iter2;
	psum += q;
	temp -= gx;
	gx = x * (a + b + iter2 - one) / (a + iter2) * gx;
	sum += q * temp;
	goto L50;
L60:
	;
    }
    ret_val = sum;

    return ret_val;
} /* ncbeta_ */

double betanc_(double x, double a, double b, double lambda,
	       double errmax, int itrmax, int *ifault)
{
    /* Initialized data */

    /* static double errmax = 1e-6; */
    /* static int itrmax = 100; */
    static double ualpha = 5.;
    static double zero = 0.;
    static double half = 0.5;
    static double one =  1.;

    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double sqrt(double), log(double), exp(double);

    /* Local variables */
    static double c__, q, a0, x0, ax, gx, xj, beta, temp, sumq, errbd;


/*     ALGORITHM AS226 APPL. STATIST. (1987) VOL. 36, NO. 2
     Incorporates modification AS R84 from AS vol. 39, pp311-2, 1990

     Returns the cumulative probability of X for the non-central beta
     distribution with parameters A, B and non-centrality LAMBDA

     Auxiliary routines required: LGAMMA - log-gamma function (ACM
     291 or AS 245), and BETAIN - incomplete-beta function (AS 63)


     Change ERRMAX and ITRMAX if desired ... */

    ret_val = x;

    *ifault = 2;
    if (lambda < zero || a <= zero || b <= zero) {
	return ret_val;
    }
    *ifault = 3;
    if (x < zero || x > one) {
	return ret_val;
    }
    *ifault = 0;
    if (x == zero || x == one) {
	return ret_val;
    }

    c__ = lambda * half;

    x0 = floor(fmax2(0., c__ - ualpha * sqrt(c__)));
    a0 = a + x0;
    beta = lgamma(a + x0) + lgamma(b) - lgamma(a0 + b);
    temp = pbeta(x, a0, b, /*lower= */TRUE, /*log_p=*/FALSE);

    gx = exp(a0 * log(x) + b * log1p( - x) - beta - log(a0));
    if (a0 > a) {
	q = exp(-c__ + x0 * log(c__)) - lgamma(x0 + 1.);
    } else {
	q = exp(-c__);
    }
    xj = zero;
    ax = q * temp;
    sumq = one - q;
    ret_val = ax;

/*     Recur over subsequent terms until convergence is achieved... */

L10:
    xj += one;
    temp -= gx;
    gx = x * (a + b + xj - one) * gx / (a + xj);
    q = q * c__ / xj;
    sumq -= q;
    ax = temp * q;
    ret_val += ax;

/*     Check for convergence and act accordingly... */

    errbd = (temp - gx) * sumq;
    if ((int) xj < itrmax && errbd > errmax) {
	goto L10;
    }
    if (errbd > errmax) {
	*ifault = 1;
    }

    return ret_val;
} /* betanc_ */

