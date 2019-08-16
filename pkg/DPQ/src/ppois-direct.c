/*
 *  DIRECT computation of ppois() {cumulative Poisson distribution function}
 *  ------
 *  instead of R's code :
 *
 *    ppois(x, lambda, low, log_p) := pgamma(lambda, x + 1, 1., !low, log_p)
 */

#ifdef DEBUG_p2
# define DEBUG_p
#endif
#ifdef DEBUG_p
# include <R_ext/Print.h>
#endif

#include "DPQpkg.h"

SEXP ppoisD(SEXP X, SEXP lambda_)
{
    if(!isReal(X))
	error("'x' must be a \"double\" numeric vector");
    double
	*x = REAL(X),
	lam = asReal(lambda_);
    R_xlen_t i, n = XLENGTH(X);
    SEXP Prob = PROTECT(allocVector(REALSXP, n));
    double *prob = REAL(Prob);

    if(ISNAN(lam)) error("lambda is NA -- invalid here");
    if(lam <= 0.) error("lambda <= 0 is invalid (== 0 : just here)");

    double
	jI = ceil(lam) - 1,
	pL = dpois(jI, lam, /*log = */FALSE);
    for(i = 0; i < n; i++) { /* prob[i] := ppois(x[i], lambda) */
	double xi = floor(x[i] + 1e-7);
	if (ISNAN(xi))
	    error("x[%d] is NA -- invalid here", i+1);
	if (xi < 0)		{ prob[i]= 0.; break; }
	if (!R_FINITE(xi))	{ prob[i]= 1.; break; }

	double f, f0;
	int j_;
	/* Compute  Prob = sum_{j = 0}^x  exp(-L) L^j / j!   {L := lambda} */
	Rboolean sml_x = (xi <= jI);
	/* the maximal term is at j_ : */
	if(sml_x) {
	    j_ = xi;
	    f = dpois(j_, lam, /*log = */FALSE);
	}
	else {
	    j_ = jI;
	    f0 = f = pL;
	}
	double P = f;
	for(int j = j_; j > 0; j--) {
	    f *= j/lam;
	    P += f;
	}
	if(!sml_x) { /* xi > ceil(lam) - 1 */
	    f = f0;
	    for(int j = j_+1; j <= xi; j++) {
		f *= lam/j;
		P += f;
	    }
	}
	prob[i] = P;
    } /* i = 1..n */
    UNPROTECT(1);
    return Prob;
}
