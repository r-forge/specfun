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

double betanc(double x, double a, double b, double lambda,
	      double errmax, int itrmax, int *ifault);

double ncbeta(double a, double b, double lambda, double x,
	      int itrmax, /* only passed to betanc() for AS 226 currently,
			   *  was hardwired = 100 */
	      double errmax,
	      int *xj, // output: gives the number of iterations taken
	      int *ifault)
{

/*       ALGORITHM AS 310 APPL. STATIST. (1997), VOL. 46, NO. 1

       Computes the cumulative distribution function of a
       non-central beta random variable
*/

    /* Builtin functions */
    double sqrt(double), log(double), exp(double), lgamma(double);
/*
      double precision function betain(x, p, q, beta, ifault)

      double precision function gammad(x, p, ifault)
*/

    // Check for admissibility of parameters :

    if (lambda <= 0. || a <= 0. || b <= 0.) {
	*ifault = 3; return x;
    }
    if (x < 0. || x > 1.) {
	*ifault = 2; return x;
    }
    if (x == 0. || x == 1.) {
	*ifault = 1; return x;
    }
    *ifault = 0;

    *xj = 0;

    if (lambda < 54.) {/* AS 226 as it stands is sufficient in this situation */

	return betanc(x, a, b, lambda, errmax, itrmax, ifault);

    } else { // lambda >= 54 :

	double q, r__, s, t, s0, s1, t0, t1, fx, gx;
	double ebd, sum, beta, temp, psum;
	int iter1, iter2;
	double errbd, ftemp;
	int iterhi, iterlo;
	double c__ = lambda * 0.5;

	int m = (int) (c__ + 0.5);
	t = 5. * sqrt((double) m);
	iterlo = m - t;
	iterhi = m + t;
	t = -c__ + m * log(c__) - lgamma(m + 1.);
	q = exp(t);
	r__ = q;
	psum = q;
	// r__1 = a + m;
	// r__2 = a + m + b;
	beta = lgamma(a + m) + lgamma(b) - lgamma(a + m + b);
	s1 = (a + m) * log(x) + b * log1p( - x) - log(a + m) - beta;
	gx = exp(s1);
	fx = gx;
	// r__1 = a + m;
/*	temp = betain_(x, &r__1, b, &beta, ifault); */
	temp = pbeta(x, a + m, b, /*lower= */TRUE, /*log_p=*/FALSE);
	ftemp = temp;
	(*xj)++;
	sum = q - temp;
	iter1 = m;

/*      The first set of iterations starts from M and goes downwards */

L20:
	if (iter1 < iterlo || q < errmax) {
	    goto L30;
	}
	q -= iter1 / c__;
	(*xj)++;
	gx = (a + iter1) / (x * (a + b + iter1 - 1.)) * gx;
	iter1--;
	temp += gx;
	psum += q;
	sum += q * temp;
	goto L20;
L30:
	t0 = lgamma(a + b) - lgamma(a + 1.) - lgamma(b);
	s0 = a * log(x) + b * log1p( - x);
	s = 0.;
	for (double j = 0.; j < iter1; ++j) {
	    s += exp(t0 + s0 + j * log(x));
	    t1 = log(a + b + j) - log(a + 1. + j) + t0;
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
	ebd = errbd + (1. - psum) * temp;
	if (ebd < errmax || iter2 >= iterhi) {
	    goto L60;
	}
	iter2 += 1.;
	(*xj)++;
	q = q * c__ / iter2;
	psum += q;
	temp -= gx;
	gx = x * (a + b + iter2 - 1.) / (a + iter2) * gx;
	sum += q * temp;
	goto L50;
L60:
	;
	return sum;
    }

} /* ncbeta */

/* call this via .C()  from R : */
void C_ncbeta(double *a, double *b, double *lambda, double *x, double *errmax,
	      int *itrmax,
	      int *ifault, double *res, int *n_res)
{
    int it_used;
    for(int i=0; i < *n_res; i++)
	res[i] = ncbeta(a[i], b[i], lambda[i], x[i],
			*itrmax, *errmax, &it_used, ifault);
}


double betanc(double x, double a, double b, double lambda,
	      double errmax, int itrmax, int *ifault)
{
/*     ALGORITHM AS226 APPL. STATIST. (1987) VOL. 36, NO. 2
     Incorporates modification AS R84 from AS vol. 39, pp311-2, 1990

     Returns the cumulative probability of X for the non-central beta
     distribution with parameters A, B and non-centrality LAMBDA

     Auxiliary routines required: LGAMMA - log-gamma function (ACM
     291 or AS 245), and BETAIN - incomplete-beta function (AS 63)
*/
    /* Builtin functions */
    double sqrt(double), log(double), exp(double);

    /* Local variables */
    static double c__, q, a0, x0, ax, gx, beta, temp, sumq, errbd;

    /* Originally

    static real errmax = 1e-6f;
    static integer itrmax = 100;

    * Change ERRMAX and ITRMAX if desired ... */

    if (lambda < 0. || a <= 0. || b <= 0.) {
	*ifault = 2; return x;
    }
    if (x < 0. || x > 1.) {
	*ifault = 3; return x;
    }
    *ifault = 0;
    if (x == 0. || x == 1.) {
	return x;
    }

    c__ = lambda * 0.5;

    x0 = floor(fmax2(0., c__ - 5. * sqrt(c__)));
    a0 = a + x0;
    beta = lgamma(a + x0) + lgamma(b) - lgamma(a0 + b);
    temp = pbeta(x, a0, b, /*lower= */TRUE, /*log_p=*/FALSE);

    gx = exp(a0 * log(x) + b * log1p( - x) - beta - log(a0));
    if (a0 > a) {
	q = exp(-c__ + x0 * log(c__)) - lgamma(x0 + 1.);
    } else {
	q = exp(-c__);
    }
    int j = 0;
    ax = q * temp;
    sumq = 1. - q;
    double ret_val = ax;

/*     Recur over subsequent terms until convergence is achieved... */
// L10:
    do {
	j += 1;
	temp -= gx;
	gx = x * (a + b + j - 1.) * gx / (a + j);
	q = q * c__ / j;
	sumq -= q;
	ax = temp * q;
	ret_val += ax;

/*     Check for convergence and act accordingly... */
	errbd = (temp - gx) * sumq;
    } while (j < itrmax && errbd > errmax);

    if (errbd > errmax) {
	*ifault = 1;
    }
    return ret_val;
} /* betanc */

