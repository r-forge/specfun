/*
 *  This is a version of R (1.9.0)'s  src/nmath/pnchisq.c
 *                                    ~~~~~~~~~~~~~~~~~~~
 *  Computes the noncentral chi-squared distribution function with
 *  positive real degrees of freedom f and nonnegative noncentrality
 *  parameter theta
 */
#include <float.h> // for DBL_*

#include "DPQpkg.h"

/* --> below:  Pnchisq_it()  called from R
 *             ------------
 * FIXME:  1) additionally vectorize !?
           2) add 'verbose' argument, instead of  '#ifdef ...' below  !!!
 */

double pnchisq_it(double x, double f, double theta,
		  double errmax, double reltol, int itrmax,
		  int *i_, int *n_terms, double *terms)
{
    double ans, lam, u, v, x2, f2, t, term, bound, f_x_2n, f_2n, lt;
    double lu = -1., l_lam = -1., l_x = -1.; /* initialized for -Wall */
    int n;
    Rboolean lamSml, tSml, is_r, is_b, is_it = FALSE;

    static const double _dbl_min_exp = M_LN2 * DBL_MIN_EXP;
    /*= -708.3964 for IEEE double precision */

    i_[0] = 0;
    *n_terms = 0;

    if (x <= 0.)	return 0.;
    if(!R_FINITE(x))	return 1.;

#ifdef DEBUG_pnch
    REprintf("pnchisq(x=%g, f=%g, theta=%g): ",x,f,theta);
#endif
    lam = .5 * theta;
    lamSml = (-lam < _dbl_min_exp);
    if(lamSml) {
	/* MATHLIB_ERROR(
	   "non centrality parameter (= %g) too large for current algorithm",
	   theta) */
        u = 0;
        lu = -lam;/* == ln(u) */
        l_lam = log(lam);
    } else {
	u = exp(-lam);
    }

    /* evaluate the first term */
    v = u;
    x2 = .5 * x;
    f2 = .5 * f;
    f_x_2n = f - x;

#ifdef DEBUG_pnch
    REprintf("-- v=exp(-th/2)=%g, x/2= %g, f/2= %g\n",v,x2,f2);
#endif

    if(f2 * DBL_EPSILON > 0.125 && /* very large f and x ~= f: probably needs */
       fabs(t = x2 - f2) <         /* other algorithm anyway */
       sqrt(DBL_EPSILON) * f2) {
	/* evade cancellation error */
	/* t = exp((1 - t)*(2 - t/(f2 + 1))) / sqrt(2*M_PI*(f2 + 1));*/
        lt = (1 - t)*(2 - t/(f2 + 1)) - 0.5 * log(2*M_PI*(f2 + 1));
#ifdef DEBUG_pnch
	REprintf(" (case I) ==> ");
#endif
    }
    else {
	/* Usual case 2: careful not to overflow .. : */
	lt = f2*log(x2) -x2 - lgammafn(f2 + 1);
    }
#ifdef DEBUG_pnch
    REprintf(" lt= %g", lt);
#endif

    tSml = (lt < _dbl_min_exp);
    if(tSml) {
	if (x > f + theta +  3* sqrt( 2*(f + 2*theta))) {
 	    /* x > E[X] + 3* sigma(X) */
	    warning("x > E[X] + 3*sigma(X) -- result may not be good");
	    // return 1.;
	    /* better than 0 --- but definitely "FIXME" */
	}
	// else {
	l_x = log(x);
	ans = term = t = 0.;
	i_[0]++; /* initial  0-term that we leave away */
	// }
    }
    else {
	t = exp(lt);
#ifdef DEBUG_pnch
 	REprintf(", t=exp(lt)= %g\n", t);
#endif
	ans = term = v * t;
    }
    terms[0] = term;

    for (n = 1, f_2n = f + 2., f_x_2n += 2.;  ; n++, f_2n += 2, f_x_2n += 2) {
#ifdef DEBUG_pnch
	REprintf("\n _OL_: n=%d",n);
#endif
	/* f_2n    === f + 2*n
	 * f_x_2n  === f - x + 2*n   > 0  <==> (f+2n)  >   x */
	if (f_x_2n > 0) {
	    /* find the error bound and check for convergence */

	    bound = t * x / f_x_2n;
#ifdef DEBUG_pnch
	    REprintf("\n L10: n=%d; term= %g; bound= %g",n,term,bound);
#endif
	    is_r = is_it = FALSE;
	    /* convergence only if BOTH absolute and relative error < 'bnd' */
	    if (((is_b = (bound <= errmax)) &&
                 (is_r = (term <= reltol * ans))) || (is_it = (n > itrmax)))
            {
#ifdef DEBUG_pnch
                REprintf("BREAK n=%d %s; bound= %g %s, rel.err= %g %s\n",
			 n, (is_it ? "> itrmax" : ""),
			 bound, (is_b ? "<= errmax" : ""),
			 term/ans, (is_r ? "<= reltol" : ""));
#endif
		break; /* out completely */
            }
	}
	else if (n > itrmax) { // (not used ?  just in case)
	    REprintf("series not converged, n=%d > itrmax (while f_x_2n = %g <= 0)\n",
		     n, f_x_2n);
	    break;
	}

	// evaluate the next term of the expansion and then the partial sum :

        if(lamSml) {
            lu += l_lam - log(n); /* u = u* lam / n */
            if(lu >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
#ifdef DEBUG_pnch
                REprintf(" n=%d; nomore underflow in u = exp(lu) ==> change\n",
			 n);
#endif
                v = u = exp(lu); /* the first non-0 'u' */
                lamSml = FALSE;
            }
        } else {
	    u *= lam / n;
	    v += u;
	}
	if(tSml) {
            lt += l_x - log(f_2n);/* t <- t * (x / f2n) */
            if(lt >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
#ifdef DEBUG_pnch
                REprintf("  n=%d; nomore underflow in t = exp(lt) ==> change\n",
			 n);
#endif
                t = exp(lt); /* the first non-0 't' */
                tSml = FALSE;
            }
        } else {
	    t *= x / f_2n;
	}
        if(!lamSml && !tSml) {
	    term = v * t;
	    ans += term;
	} else {
	    i_[0]++; /* another initial 0-term */
	}
	terms[n] = term;

    } /* for(n ...) */

    if (is_it) {
	REprintf("pnchisq(x=%g, ..): not converged in %d iter.",
		 x, itrmax);
    }
#ifdef DEBUG_pnch
    REprintf("\n == L_End: n=%d; term= %g; bound=%g\n",n,term,bound);
#endif
    *n_terms = n;
    return (ans);
}

// see 'FIXME' above !!
void Pnchisq_it(double *x, double *f, double *theta,
		double *errmax, double *reltol, int *itrmax,
		int *i_0, int *n_terms, double *terms, double *prob)
{
    *prob = pnchisq_it(*x, *f, *theta, *errmax, *reltol, *itrmax,
		       i_0, n_terms, terms);
}
