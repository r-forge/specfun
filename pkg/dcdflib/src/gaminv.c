/* The NSWC library (1993) has an improved version of gaminv()
 *     ------------
 * --> /usr/local/lib.src/NSWC/nswc/
 * and /usr/local/lib.src/NSWC/C-src/ erfi-1.c, pni-1.c
 *
 * I have converted them via f2c & f2c-clean and by hand:
 *			Martin Maechler, Statistik, ETH Zurich, May 1999
 *
 * The `older' version in dcdflib (1.1; Dec.1997) is now in
 * ./dcdf-gaminv.c
 */

/* Produced by
 * $Id: f2c-clean,v 1.6 1999/05/17 17:11:30 maechler Exp $

 * erfi.f, pni.f, gaminv.f -- translated by f2c (version 19980913).
 */

#include <R.h>
#include <Rmath.h>

#include "cdflib.h"

/* Table of constant values */

static int I0 = 0;

void gaminv(double *a, double *x, double *x0,
	    double *p, double *q, int *ierr)
{
/*
 ----------------------------------------------------------------------
	    INVerse incomplete GAMma ratio function

	Given positive a, and nonegative p and q where p + q = 1,
	then x is computed where p(a,x) = p and q(a,x) = q.
	Schroder iteration is employed.
	The routine attempts to compute x to 10 significant digits
					  ~~~~~~~~~~~~~~~~~~~~~~~~
	if this is possible for the particular computer arithmetic being used.

			 ------------

	x is a variable. if p = 0 then x is assigned the value 0,
	and if q = 0 then x is set to the largest floating point
	number available.  Otherwise, gaminv attempts to obtain
	a solution for p(a,x) = p and q(a,x) = q.
	If the routine is successful then the solution is stored in x.

	x0:  optional initial approximation for x.
	     if the user does not wish to supply an initial approximation,
	     then set x0 <= 0.

	ierr:   a variable that reports the status of the results.
	        when the routine terminates, ierr has one of the following
	        values :

	   ierr =  0	the solution was obtained. iteration was not used.
	   ierr >  0	the solution was obtained. ierr iterations
			were performed.
	   ierr = -2	(input error) a <= 0
	   ierr = -3	no solution was obtained. the ratio q/a
			is too large.
	   ierr = -4	(input error) p or q is negative, or  p + q != 1.
	   ierr = -6	20 iterations were performed. the most
			recent value obtained for x is given.
			this cannot occur if x0 <= 0.
	   ierr = -7	iteration failed. no value is given for x.
			this may occur when x is approximately 0.
	   ierr = -8	a value for x has been obtained, but the
			routine is not certain of its accuracy.
			iteration cannot be performed in this
			case. if x0 <= 0, this can occur only
			when p or q is approximately 0. if x0 is
			positive then this can occur when a is
			exceedingly close to x and a is extremely
			large (say a >= 1.e20).
 ----------------------------------------------------------------------
     written by Alfred H. Morris, Jr.
	Naval Surface Warfare Center (NSWC)
	Dahlgren, Virginia
     -------------------
     Revised ... January 1992
*/

    /*static double ln10 = 2.302585; -- use M_LN10 */ /* ln10 */
    static double cE = .577215664901533;/* = Euler's Constant */
    static double bmin[2] = { 1e-28,1e-13 };
    static double emin[2] = { .002, .006 };
    static double tol = 1e-5;

    /* System generated locals */
    double r__1;

    /* Local variables */
    double amin, amax, xmin, eps, rta, sum, b = 0./*Wall*/;
    double d__, e, g, h__, s, t, u, r__, w, y, z__, c1, c2, c3, e2, c4, c5;
    double s2, qg, pn, qn, xn, am1,ap1, ap2, ap3, apn;
    long ier, iop;

/*
  Machine dependent constants :
	e is the smallest number for which 1 + e > 1,
	xmin is the smallest positive number.
*/
    e    = spmpar(1);
    xmin = spmpar(2);

    /* Argument checking: */
    if (*a <= 0.) { *ierr = -2; return; }
    if (*p < 0. || *q < 0.) { *ierr = -4; return; }
    t = *p + *q - .5 - .5;
    if (fabs(t) > fmax2(e,1e-15) * 5.) { *ierr = -4; return; }

    *x = 0.;
    *ierr = 0;
    xmin /= e;
    if (*p / e <= xmin)	goto L400;
    if (*q / e <= xmin) goto L560;

    if (*a == 1.)	goto L410;

    e2 = e + e;
    amax = 4e-11 / (e * e);
    eps = fmax2(100 * e, 1e-10);

    iop = 1;
    if (e > 1e-10) {
	iop = 2;
    }
    xn = *x0;
    if (*x0 > 0.)	goto L100;

/*        SELECTION OF THE INITIAL APPROXIMATION XN OF X
                       WHEN A < 1 */

    if (*a > 1.) {
	goto L50;
    }
    r__1 = *a + 1.;
    g = Xgamm(&r__1);
    qg = *q * g;
    if (qg == 0.) {
	goto L560;
    }
    b = qg / *a;
    if (qg > *a * .6) {
	goto L20;
    }
    if (*a >= .3 || b < .35) {
	goto L10;
    }
    t = exp(-(b + cE));
    u = t * exp(t);
    xn = t * exp(u);
    goto L100;

L10:
    if (b >= .45) {
	goto L20;
    }
    if (b == 0.) {
	goto L560;
    }
    y = -log(b);
    s = .5 - *a + .5;
    z__ = log(y);
    t = y - s * z__;
    if (b >= .15) {
	xn = y - s * log(t) - log(s / (t + 1.) + 1.);
	goto L200;
    }
/*L11:*/
    if (b > .01) {
	u = ((t + (3. - *a) * 2.) * t + (2. - *a) * (3. - *a)) /
	    ((t + (5. - *a)) * t + 2.);
	xn = y - s * log(t) - log(u);
	goto L200;
    }
Loop12:
    c1 = -s * z__;
    c2 = -s * (c1 + 1.);
    c3 = s * ((c1 * .5 + (2. - *a)) * c1 + (2.5 - *a * 1.5));
    c4 = -s * (((c1 / 3. + (2.5 - *a * 1.5)) * c1 + ((*a - 6.) * *a + 7.)
	    ) * c1 + ((*a * 11. - 46.) * *a + 47.) / 6.);
    c5 = -s * ((((-c1 / 4. + (*a * 11. - 17.) / 6.) * c1 + ((*a * -3. +

	    13.) * *a - 13.)) * c1 + (((*a * 2. - 25.) * *a + 72.) * *a

	    - 61.) * .5) * c1 + (((*a * 25. - 195.) * *a + 477.) * *a -

	    379.) / 12.);
    xn = (((c5 / y + c4) / y + c3) / y + c2) / y + c1 + y;
    if (*a > 1.) {
	goto L200;
    }
    if (b > bmin[iop - 1]) {
	goto L200;
    }
    *x = xn;
    return;

L20:
    if (b * *q > 1e-8) {
	goto L21;
    }
    xn = exp(-(*q / *a + cE));
    goto L30;
L21:
    if (*p <= .9) {
	goto L22;
    }
    r__1 = -(*q);
    xn = exp((alnrel(&r__1) + gamln1(a)) / *a);
    goto L30;
L22:
    xn = exp(log(*p * g) / *a);

L30:
    if (xn == 0.) {
	goto L510;
    }
    t = .5 - xn / (*a + 1.) + .5;
    xn /= t;
    goto L100;

/*        SELECTION OF THE INITIAL APPROXIMATION XN OF X
                       WHEN A > 1 */

L50:
    t = *p - .5;
    if (*q < .5) {
	t = .5 - *q;
    }
    pni(p, q, &t, &s, &ier);

    rta = sqrt(*a);
    s2 = s * s;
    xn = (((s2 * 12. - 243.) * s2 - 923.) * s2 + 1472.) / 204120.;
    xn = xn / *a + s * ((s2 * 9. + 256.) * s2 - 433.) / (rta * 38880.) - (
	    (s2 * 3. + 7.) * s2 - 16.) / 810.;
    xn = *a + s * rta + (s2 - 1.) / 3. + s * (s2 - 7.) / (rta * 36.) + xn

	    / *a;
    xn = fmax2(xn,0.);

    amin = 20.;
    if (e < 1e-8) {
	amin = 250.;
    }
    if (*a < amin) {
	goto L60;
    }
    *x = xn;
    d__ = .5 - *x / *a + .5;
    if (fabs(d__) <= .1) {
	return;
    }

L60:
    if (*p <= .5) {
	goto L70;
    }
    if (xn < *a * 3.) {
	goto L200;
    }
    w = log(*q);
    y = -(w + gamln(a));
   d__ = fmax2(2, *a * (*a - 1.));
    if (y < M_LN10 * d__)
	goto L61;

    s = 1. - *a;
    z__ = log(y);
    goto Loop12;

L61:
    t = *a - 1.;
    r__1 = -t / (xn + 1.);
    xn = y + t * log(xn) - alnrel(&r__1);
    goto L200;

L70:
    ap1 = *a + 1.;
    if (xn > ap1 * .7) {
	goto L101;
    }
    w = log(*p) + gamln(&ap1);
    if (xn > ap1 * .15) {
	goto L80;
    }
    ap2 = *a + 2.;
    ap3 = *a + 3.;
    *x = exp((w + *x) / *a);
    *x = exp((w + *x - log(*x / ap1 * (*x / ap2 + 1.) + 1.)) / *a);
    *x = exp((w + *x - log(*x / ap1 * (*x / ap2 + 1.) + 1.)) / *a);
    *x = exp((w + *x - log(*x / ap1 * (*x / ap2 * (*x / ap3 + 1.) + 1.) +

	    1.)) / *a);
    xn = *x;
    if (xn > ap1 * .01) {
	goto L80;
    }
    if (xn <= emin[iop - 1] * ap1) {
	return;
    }
    goto L101;

L80:
    apn = ap1;
    t = xn / apn;
    sum = t + 1.;
/*L81:*/
    do {
	apn += 1.;
	t *= xn / apn;
	sum += t;
    } while (t > 1e-4);

    t = w - log(sum);
    xn = exp((xn + t) / *a);
    xn *= 1. - (*a * log(xn) - xn - t) / (*a - xn);
    goto L101;

/*                 SCHRODER ITERATION USING P */

L100:
    if (*p > .5) {
	goto L200;
    }
L101:
    if (*p <= xmin) {
	goto L550;
    }
    am1 = *a - .5 - .5;
L102:
    if (*a <= amax) {
	goto L110;
    }
    d__ = .5 - xn / *a + .5;
    if (fabs(d__) <= e2) {
	goto L550;
    }

L110:
    if (*ierr >= 20) {
	goto L530;
    }
    ++(*ierr);
    gratio(a, &xn, &pn, &qn, &I0);
    if (pn == 0. || qn == 0.) {
	goto L550;
    }
    r__ = rcomp(a, &xn);
    if (r__ < xmin) {
	goto L550;
    }
    t = (pn - *p) / r__;
    w = (am1 - xn) * .5;
    if (fabs(t) <= .1 && fabs(w * t) <= .1) {
	goto L120;
    }
    *x = xn * (1. - t);
    if (*x <= 0.) {
	goto L540;
    }
    d__ = fabs(t);
    goto L121;

L120:
    h__ = t * (w * t + 1.);
    *x = xn * (1. - h__);
    if (*x <= 0.) {
	goto L540;
    }
    if (fabs(w) >= 1. && fabs(w) * t * t <= eps) {
	return;
    }
    d__ = fabs(h__);
L121:
    xn = *x;
    if (d__ > tol) {
	goto L102;
    }
    if (d__ <= eps) {
	return;
    }
    if (fabs(*p - pn) <= tol * *p) {
	return;
    }
    goto L102;

/*                 SCHRODER ITERATION USING Q */

L200:
    if (*q <= xmin) {
	goto L550;
    }
    am1 = *a - .5 - .5;
L201:
    if (*a <= amax) {
	goto L210;
    }
    d__ = .5 - xn / *a + .5;
    if (fabs(d__) <= e2) {
	goto L550;
    }

L210:
    if (*ierr >= 20) {
	goto L530;
    }
    ++(*ierr);
    gratio(a, &xn, &pn, &qn, &I0);
    if (pn == 0. || qn == 0.) {
	goto L550;
    }
    r__ = rcomp(a, &xn);
    if (r__ < xmin) {
	goto L550;
    }
    t = (*q - qn) / r__;
    w = (am1 - xn) * .5;
    if (fabs(t) <= .1 && fabs(w * t) <= .1) {
	goto L220;
    }
    *x = xn * (1. - t);
    if (*x <= 0.) {
	goto L540;
    }
    d__ = fabs(t);
    goto L221;

L220:
    h__ = t * (w * t + 1.);
    *x = xn * (1. - h__);
    if (*x <= 0.) {
	goto L540;
    }
    if (fabs(w) >= 1. && fabs(w) * t * t <= eps) {
	return;
    }
    d__ = fabs(h__);
L221:
    xn = *x;
    if (d__ > tol) {
	goto L201;
    }
    if (d__ <= eps) {
	return;
    }
    if (fabs(*q - qn) <= tol * *q) {
	return;
    }
    goto L201;

/*                       SPECIAL CASES */

L400:
    *ierr = -8;
    return;

L410:
    if (*q < .9) {
	goto L411;
    }
    r__1 = -(*p);
    *x = -alnrel(&r__1);
    return;
L411:
    *x = -log(*q);
    return;

/*                       ERROR RETURN */



L510:    *ierr = -3; return;

/*L520:    *ierr = -4; return;*/

L530:    *ierr = -6; return;

L540:    *ierr = -7; return;

L550:    *x = xn; *ierr = -8; return;

L560:    *x = spmpar(3); *ierr = -8; return;
} /* gaminv() */


double erfi(double *p, double *q)
{
    /* Initialized data */

    static double c__ = .5625;
    static double a3[9] = { 3.20540542206205e-9,1.899479322632128e-6,
	    2.814223189858532e-4,.01370504879067817,.2268143542005976,
	    1.09842195989234,.6791143397056208,-.834334189167721,
	    .3421951267240343 };
    static double b3[6] = { 3.205405053282398e-9,1.899480592260143e-6,
	    2.81434969109894e-4,.01371092249602266,.2275172815174473,
	    1.125348514036959 };
    static double c1 = .87890625;
    static double c2 = -230.2585092994046;
    static double a[6] = { 140.0216916161353,-720.4275515686407,
	    1296.708621660511,-969.7932901514031,276.2427049269425,
	    -20.12940180552054 };
    static double b[6] = { 129.1046303114685,-731.2308064260973,
	    1494.970492915789,-1337.793793683419,503.3747142783567,
	    -62.20205554529216 };
    static double a1[7] = { -.1690478046781745,3.524374318100228,
	    -26.98143370550352,93.40783041018743,-145.5364428646732,
	    88.05852004723659,-13.49018591231947 };
    static double b1[7] = { -.1203221171313429,2.684812231556632,
	    -22.42485268704865,87.23495028643494,-160.4352408444319,
	    125.9117982101525,-31.84861786248824 };
    static double a2[9] = { 3.100808562552958e-5,.00409748760301194,
	    .1214902662897276,1.109167694639028,3.228379855663924,
	    2.881691815651599,2.047972087262996,.8545922081972148,
	    .003551095884622383 };
    static double b2[8] = { 3.100809298564522e-5,.004097528678663915,
	    .1215907800748757,1.118627167631696,3.43236398430529,
	    4.140284677116202,4.119797271272204,2.162961962641435 };


    /* Local variables */
    double s, t, v, v1, eps;

/* -----------------------------------------------------------------------

              EVALUATION OF THE INVERSE ERROR FUNCTION

                        ---------------

     FOR 0 <= P < 1,  W = ERFI(P,Q) WHERE ERF(W) = P. IT IS
     ASSUMED THAT Q = 1 - P.  IF P < 0, Q <= 0, OR P + Q IS
     NOT 1, THEN ERFI(P,Q) IS SET TO A NEGATIVE VALUE.

 -----------------------------------------------------------------------
     REFERENCE. MATHEMATICS OF COMPUTATION,OCT.1976,PP.827-830.
                  J.M.BLAIR,C.A.EDWARDS,J.H.JOHNSON
 -----------------------------------------------------------------------
 -----------------------------------------------------------------------
     C2 = LN(1.E-100)
 -----------------------------------------------------------------------
 -----------------------------------------------------------------------
                            TABLE NO.16
 -----------------------------------------------------------------------
 -----------------------------------------------------------------------
                            TABLE NO.36
 -----------------------------------------------------------------------
 -----------------------------------------------------------------------
                            TABLE NO.56
 -----------------------------------------------------------------------
 -----------------------------------------------------------------------
                            TABLE NO.79
 -----------------------------------------------------------------------
 ----------------------------------------------------------------------- */
    if (*p < 0. || *q <= 0.) {
	goto L100;
    }

    eps = fmax2(spmpar(1),1e-15);
    t = .5 - (*p + *q) + .5;
    if (fabs(t) > eps * 3.) {
	goto L110;
    }

/*                      0 <= P <= 0.75 */

    if (*p > .75) {
	goto L10;
    }
    v = *p * *p - c__;
    t = *p * (((((a[5] * v + a[4]) * v + a[3]) * v + a[2]) * v + a[1]) * v +

	    a[0]);
    s = (((((v + b[5]) * v + b[4]) * v + b[3]) * v + b[2]) * v + b[1]) * v +

	    b[0];
    goto L40;

/*                    0.75 < P <= 0.9375 */

L10:
    if (*p > .9375) {
	goto L20;
    }
    v = *p * *p - c1;
    t = *p * ((((((a1[6] * v + a1[5]) * v + a1[4]) * v + a1[3]) * v + a1[2]) *
	     v + a1[1]) * v + a1[0]);
    s = ((((((v + b1[6]) * v + b1[5]) * v + b1[4]) * v + b1[3]) * v + b1[2]) *
	     v + b1[1]) * v + b1[0];
    goto L40;

/*                  1.E-100 <= Q < 0.0625 */

L20:
    v1 = log(*q);
    v = 1. / sqrt(-v1);
    if (v1 < c2) {
	goto L30;
    }
    t = (((((((a2[8] * v + a2[7]) * v + a2[6]) * v + a2[5]) * v + a2[4]) * v

	    + a2[3]) * v + a2[2]) * v + a2[1]) * v + a2[0];
    s = v * ((((((((v + b2[7]) * v + b2[6]) * v + b2[5]) * v + b2[4]) * v +

	    b2[3]) * v + b2[2]) * v + b2[1]) * v + b2[0]);
    goto L40;

/*                 1.E-10000 <= Q < 1.E-100 */

L30:
    t = (((((((a3[8] * v + a3[7]) * v + a3[6]) * v + a3[5]) * v + a3[4]) * v

	    + a3[3]) * v + a3[2]) * v + a3[1]) * v + a3[0];
    s = v * ((((((v + b3[5]) * v + b3[4]) * v + b3[3]) * v + b3[2]) * v + b3[
	    1]) * v + b3[0]);
L40:
    return  t / s;

/*                         ERROR RETURN */

L100:    return -1.;

L110:    return -2.;

} /* erfi */

void pni(double *p, double *q, double *d, double *w, long *ierr)
{
/* -----------------------------------------------------------------------

         EVALUATION OF THE INVERSE NORMAL DISTRIBUTION FUNCTION

                           ------------

     LET F(T) = 1/(SQRT(2*PI)*EXP(-T*T/2)). THEN THE FUNCTION

        PROB(X) = INTEGRAL FROM MINUS INFINITY TO X OF F(T)

     IS THE NORMAL DISTRIBUTION FUNCTION OF ZERO MEAN AND UNIT
     VARIANCE. IT IS ASSUMED THAT P > 0, Q > 0, P + Q = 1,
     AND D = P - 0.5. THE VALUE W IS COMPUTED WHERE PROB(W) = P.

     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.

       IERR = 0  NO INPUT ERRORS WERE DETECTED. W WAS COMPUTED.
       IERR = 1  EITHER P OR Q IS INCORRECT.
       IERR = 2  D IS INCORRECT.

 -----------------------------------------------------------------------
*/

    /*static double rt2 = 1.414213562373095; -- use M_SQRT2  = SQRT(2) */
    double t, u, v, eps;

    t = fmin2(*p,*q);
    if (t <= 0.) {
	goto Err_1;
    }

    eps = fmax2(spmpar(1),1e-15);
    *w = .5 - (*p + *q) + .5;
    if (fabs(*w) > eps * 2.) {
	goto Err_1;
    }

    u = fabs(*d + *d);
    v = t + t;
    *w = erfi(&u, &v);
    if (*w < 0) {
	*ierr = 2; return;
    }

    *w = M_SQRT2 * *w;
    if (*d < 0)	*w = -(*w);
    *ierr = 0;    return;

Err_1: *ierr = 1; return;

} /* pni() */

