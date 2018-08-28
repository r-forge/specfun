/*
 *  DIRECT computation of ppois() {cumulative Poisson distribution function}
 *  ------
 *  instead of R's code :
 *
 *    ppois(x, lambda, low, log_p) := pgamma(lambda, x + 1, 1., !low, log_p)
 */

#include <math.h>

#include <R_ext/Arith.h>
#include <R_ext/Error.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

#ifdef DEBUG_p2
# define DEBUG_p
#endif
#ifdef DEBUG_p
# include <R_ext/Print.h>
#endif

#include "DPQpkg.h"

void ppois_D(double *X, int *n, double *lambda, /* result: */ double *prob)
{
    double x, lam, f, pL, f0, P;
    Rboolean sml_x;
    int i, j, jI, j_;

    if (ISNAN(*lambda))
	error("lambda is NA -- invalid here");

    lam = *lambda;
    if(lam <= 0.) error("lambda <= 0 is invalid (== 0 : just here");

    jI = ceil(lam) - 1;
    pL = dpois(jI, lam, /*log = */FALSE);
    for(i = 0; i < *n; i++) { /* prob[i] := ppois(x[i], lambda) */
	x = floor(X[i] + 1e-7);
	if (ISNAN(x))
	    error("x[%d] is NA -- invalid here", i+1);
	if (x < 0)		{ prob[i]= 0.; break; }
	if (!R_FINITE(x))	{ prob[i]= 1.; break; }

	/* Compute  Prob = sum_{j = 0}^x  exp(-L) L^j / j!   {L := lambda} */
	sml_x = (x <= jI);
	/* the maximal term is at j_ : */
	if(sml_x) {
	    j_ = x;
	    f = dpois(j_, lam, /*log = */FALSE);
	}
	else {
	    j_ = jI;
	    f0 = f = pL;
	}
	P = f;
	for(j = j_; j > 0; j--) {
	    f *= j/lam;
	    P += f;
	}
	if(!sml_x) { /* x > ceil(lam) - 1 */
	    f = f0;
	    for(j = j_+1; j <= x; j++) {
		f *= lam/j;
		P += f;
	    }
	}
	prob[i] = P;
    } /* i = 1..n */
    return;
}
