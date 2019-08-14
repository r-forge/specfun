/* 310-CumNonCentralBeta.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

doublereal ncbeta_(real *a, real *b, real *lambda, real *x, real *errmax, 
	integer *ifault)
{
    /* Initialized data */

    static real zero = 0.f;
    static real half = .5f;
    static real one = 1.f;
    static real five = 5.f;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static real c__;
    static integer i__, j, m;
    static real q, r__, s, t, s0, s1, t0, t1, fx, gx;
    static integer xj;
    static real ebd, sum, beta, temp, psum;
    static integer iter1, iter2;
    static real errbd, ftemp;
    extern doublereal gammad_(real *, real *, integer *), betanc_(real *, 
	    real *, real *, real *, integer *), alngam_(real *), betain_(real 
	    *, real *, real *, real *, integer *), double_(integer *);
    static integer iterhi, iterlo;


/*       ALGORITHM AS 310 APPL. STATIST. (1997), VOL. 46, NO. 1 */

/*       Computes the cumulative distribution function of a */
/*       non-central beta random variable */



/*       Local variable XJ gives the number of iterations taken */




    ret_val = *x;

/*       Check for admissibility of parameters */

    *ifault = 3;
    if (*lambda <= zero || *a <= zero || *b <= zero) {
	return ret_val;
    }
    *ifault = 2;
    if (*x < zero || *x > one) {
	return ret_val;
    }
    *ifault = 1;
    if (*x == zero || *x == one) {
	return ret_val;
    }
    *ifault = 0;

    c__ = *lambda * half;
    xj = zero;

    if (*lambda < 54.f) {

/*       AS 226 as it stands is sufficient in this situation */

	ret_val = betanc_(x, a, b, lambda, ifault);
	return ret_val;
    } else {
	m = (integer) (c__ + half);
	t = five * sqrt(double_(&m));
	iterlo = m - t;
	iterhi = m + t;
	r__1 = m + one;
	t = -c__ + m * log(c__) - alngam_(&r__1);
	q = exp(t);
	r__ = q;
	psum = q;
	r__1 = *a + m;
	r__2 = *a + m + *b;
	beta = alngam_(&r__1) + alngam_(b) - alngam_(&r__2);
	s1 = (*a + m) * log(*x) + *b * log(one - *x) - log(*a + m) - beta;
	gx = exp(s1);
	fx = gx;
	r__1 = *a + m;
	temp = betain_(x, &r__1, b, &beta, ifault);
	ftemp = temp;
	xj += one;
	sum = q - temp;
	iter1 = m;

/*      The first set of iterations starts from M and goes downwards */

L20:
	if (iter1 < iterlo) {
	    goto L30;
	}
	if (q < *errmax) {
	    goto L30;
	}
	q -= iter1 / c__;
	xj += one;
	gx = (*a + iter1) / (*x * (*a + *b + iter1 - one)) * gx;
	iter1 -= one;
	temp += gx;
	psum += q;
	sum += q * temp;
	goto L20;
L30:
	r__1 = *a + *b;
	r__2 = *a + one;
	t0 = alngam_(&r__1) - alngam_(&r__2) - alngam_(b);
	s0 = *a * log(*x) + *b * log(one - *x);

	i__1 = iter1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = i__ - one;
	    s += exp(t0 + s0 + j * log(*x));
	    t1 = log(*a + *b + j) - log(*a + one + j) + t0;
	    t0 = t1;
/* L40: */
	}

/*       Compute the first part of error bound */

	r__1 = (real) iter1;
	errbd = (one - gammad_(&c__, &r__1, ifault)) * (temp + s);
	q = r__;
	temp = ftemp;
	gx = fx;
	iter2 = m;
L50:
	ebd = errbd + (one - psum) * temp;
	if (ebd < *errmax || iter2 >= iterhi) {
	    goto L60;
	}
	iter2 += one;
	xj += one;
	q = q * c__ / iter2;
	psum += q;
	temp -= gx;
	gx = *x * (*a + *b + iter2 - one) / (*a + iter2) * gx;
	sum += q * temp;
	goto L50;
L60:
	;
    }
/* L70: */
    ret_val = sum;

    return ret_val;
} /* ncbeta_ */

doublereal betanc_(real *x, real *a, real *b, real *lambda, integer *ifault)
{
    /* Initialized data */

    static real errmax = 1e-6f;
    static integer itrmax = 100;
    static real ualpha = 5.f;
    static real zero = 0.f;
    static real half = .5f;
    static real one = 1.f;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static real c__, q, a0, x0, ax, gx, xj, beta, temp, sumq, errbd;
    extern doublereal alogam_(real *, integer *), betain_(real *, real *, 
	    real *, real *, integer *);


/*     ALGORITHM AS226 APPL. STATIST. (1987) VOL. 36, NO. 2 */
/*     Incorporates modification AS R84 from AS vol. 39, pp311-2, 1990 */

/*     Returns the cumulative probability of X for the non-central beta */
/*     distribution with parameters A, B and non-centrality LAMBDA */

/*     Auxiliary routines required: ALOGAM - log-gamma function (ACM */
/*     291 or AS 245), and BETAIN - incomplete-beta function (AS 63) */


/*     Change ERRMAX and ITRMAX if desired ... */



    ret_val = *x;

    *ifault = 2;
    if (*lambda < zero || *a <= zero || *b <= zero) {
	return ret_val;
    }
    *ifault = 3;
    if (*x < zero || *x > one) {
	return ret_val;
    }
    *ifault = 0;
    if (*x == zero || *x == one) {
	return ret_val;
    }

    c__ = *lambda * half;

/*     Initialize the series ... */

/* Computing MAX */
    r__1 = c__ - ualpha * sqrt(c__);
    x0 = (real) ((integer) dmax(r__1,zero));
    a0 = *a + x0;
    r__1 = a0 + *b;
    beta = alogam_(&a0, ifault) + alogam_(b, ifault) - alogam_(&r__1, ifault);
    temp = betain_(x, &a0, b, &beta, ifault);
    gx = exp(a0 * log(*x) + *b * log(one - *x) - beta - log(a0));
    if (a0 > *a) {
	r__1 = x0 + one;
	q = exp(-c__ + x0 * log(c__)) - alogam_(&r__1, ifault);
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
    gx = *x * (*a + *b + xj - one) * gx / (*a + xj);
    q = q * c__ / xj;
    sumq -= q;
    ax = temp * q;
    ret_val += ax;

/*     Check for convergence and act accordingly... */

    errbd = (temp - gx) * sumq;
    if ((integer) xj < itrmax && errbd > errmax) {
	goto L10;
    }
    if (errbd > errmax) {
	*ifault = 1;
    }

    return ret_val;
} /* betanc_ */

