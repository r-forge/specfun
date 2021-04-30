/*
 *  Miscellanous  C  functions for pkg  'DPQ'
 *  ------
 */

// #include <float.h> /* DBL_MIN etc */

// for error() or if(verbose) ..
#include <R_ext/Print.h>

#include "DPQpkg.h"


// To be  .Call()ed  from R : ---------------

// Computes 'log(1 + X) - X'  accurately even for  |x| << 1
SEXP R_log1pmx(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ?
				Rf_duplicate(x_) : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n);
    double *x = REAL(x_), *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = log1pmx(x[i]); // log1pmx() implemented in R's src/nmath/pgamma.c
    UNPROTECT(1);
    return r_;
}


// Computes 'log(1 + exp(X))' accurately, notably for large x, e.g., x > 720.
SEXP R_log1pexp(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ?
				Rf_duplicate(x_) : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n);
    double *x = REAL(x_), *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = log1pexp(x[i]);
    UNPROTECT(1);
    return r_;
}

#include <Rversion.h>
#if R_VERSION < R_Version(4,1,0)
// From --- src/nmath/{dpq.h, plogis.c} ---- :
// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
double log1mexp(double x) { return R_Log1_Exp(-x); }
#endif


/** Computes 'log(1 - exp(-x))' accurately, carefully for two regions of x,
 * optimally cutting off at log 2 (= 0.693147..),
 * if(x <= log(2)  log(-expm1(-x))  else  log1p(-exp(-x))
 */
SEXP R_log1mexp(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ?
				Rf_duplicate(x_) : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n);
    double *x = REAL(x_), *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = log1mexp(x[i]);
    UNPROTECT(1);
    return r_;
}


// Computes 'log(gamma(X + 1))'  accurately even for small x, i.e., 0 < x < 0.5.
SEXP R_lgamma1p(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ?
				Rf_duplicate(x_) : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n);
    double *x = REAL(x_), *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = lgamma1p(x[i]);
    UNPROTECT(1);
    return r_;
}

/* in C have     r = frexp (x, &e);
 * here, as "proper function" of _1_ argument, returning "two results" as list */
SEXP R_frexp(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ?
				 Rf_duplicate(x_) : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n), // (r_, e_) will be protected as parts of ans
	 e_ = allocVector(INTSXP,  n),
	ans = PROTECT(allocVector(VECSXP, 2)), nms; // --> list(r = ., e = .)
    SET_VECTOR_ELT(ans, 0, r_);
    SET_VECTOR_ELT(ans, 1, e_);
    setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, mkChar("r"));
    SET_STRING_ELT(nms, 1, mkChar("e"));

    double *x = REAL(x_), *r = REAL(r_);
    int *e = INTEGER(e_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = frexp(x[i], e+i);
    UNPROTECT(2);
    return ans;
}

// ldexp(f, E) := f * 2^E
SEXP R_ldexp(SEXP f_, SEXP E_)
{
    R_xlen_t n = XLENGTH(PROTECT(f_ = isReal(f_) ?
				 Rf_duplicate(f_) : coerceVector(f_, REALSXP)));
    PROTECT(E_ = isInteger(E_) ? Rf_duplicate(E_) : coerceVector(E_, INTSXP));
    if(XLENGTH(E_) != n) // FIXME(?) --> recycling (f, E) !
	error(_("'E' is not of the same length as 'f': %.0f != %.0f"), (double)n, (double)XLENGTH(E_));
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    double *f = REAL(f_), *r = REAL(r_);
    int *E = INTEGER(E_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = ldexp(f[i], E[i]);
    UNPROTECT(3);
    return r_;
}
