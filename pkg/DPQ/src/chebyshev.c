/* Apart from these first lines: DIRECT COPY  from R sources  of
   <R>/src/nmath/chebyshev.c
       --------------------- @MM = ~/R/D/r-devel/R/src/nmath/chebyshev.c  */
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    int chebyshev_init(double *dos, int nos, double eta)
 *    double chebyshev_eval(double x, double *a, int n)
 *
 *  DESCRIPTION
 *
 *    "chebyshev_init" determines the number of terms for the
 *    double precision orthogonal series "dos" needed to insure
 *    the error is no larger than "eta".  Ordinarily eta will be
 *    chosen to be one-tenth machine precision.
 *
 *    "chebyshev_eval" evaluates the n-term Chebyshev series
 *    "a" at "x".
 *
 *  NOTES
 *
 *    These routines are translations into C of Fortran routines
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    Based on the Fortran routine dcsevl by W. Fullerton.
 *    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
 */

#include "DPQpkg.h"

/* NaNs propagated correctly */


int chebyshev_init(const double dos[], int nos, double eta)
{
    if (nos < 1)
	return 0;

    double err = 0.0;
    int i = 0;			/* just to avoid compiler warnings */
    for (int ii=1; ii <= nos; ii++) {
	i = nos - ii;
	err += fabs(dos[i]);
	if (err > eta) {
	    return i;
	}
    }
    return i;
}


// Chebyshev Polynomial P(x; a[]), x in [-1,1], coefficients a[0..(n-1)]
double chebyshev_eval(double x, const double a[], const int n)
{
/* In 'DPQ' pkg:  disable these to allow for more experimentation, incl. *extra*polation
    if (n < 1 || n > 1000) ML_WARN_return_NAN;

    if (x < -1.1 || x > 1.1) ML_WARN_return_NAN;
*/
    double
	twox = x * 2,
	b2 = 0, b1 = 0, b0 = 0;
    for (int i = 1; i <= n; i++) { // i in 1:n  <==>  n-i  in  0:(n-1)
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

// To be  .Call()ed  from R :
SEXP R_chebyshev_eval(SEXP x_, SEXP a_, SEXP n_)
{
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    PROTECT(a_ = isReal(a_) ? a_ : coerceVector(a_, REALSXP));
    // n_ is checked via asInteger() below
    R_xlen_t i,	nx = XLENGTH(x_);
    SEXP r_ = PROTECT(allocVector(REALSXP, nx));
    double *x = REAL(x_), *a = REAL(a_), *r = REAL(r_);
    int ntrm = asInteger(n_);
    if(ntrm <= 0)
	error("ntrm = %d <= 0", ntrm);
    if(ntrm > LENGTH(a_))
	error("ntrm = %d > length(a) = %d", ntrm, LENGTH(a_));
    for(i=0; i < nx; i++) {
	r[i] = chebyshev_eval(x[i], a, ntrm);
    }
    UNPROTECT(3);
    return r_;
}

SEXP R_chebyshev_init(SEXP coef_, SEXP eta_)
{
    PROTECT(coef_ = isReal(coef_) ? coef_ : coerceVector(coef_, REALSXP));
    if(XLENGTH(coef_) > INT_MAX)
	error("length(coef) = %ld > max.int = %d", (long)XLENGTH(coef_), INT_MAX);
    int ntrms = chebyshev_init(REAL(coef_), LENGTH(coef_), asReal(eta_));
    UNPROTECT(1);
    return ScalarInteger(ntrms);
}
