#include "DPQpkg.h"

/*---------------------------------------------vvvvvvvvvvvvvvvvvvvvv-------
 * Cut'n'paste from R's sources  ~/R/D/r-devel/R/src/nmath/pgamma.c
 * [as of 2021-04-19, lines 67--118]
 * ADDING 'maxit' trace option
 *-------------------------------------------------------------------------*/


/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxiliary in log1pmx() and lgamma1p()
 */
static double
logcf (double x, double i, double d, double eps /* ~ relative tolerance */,
       int trace)
{
    const int maxit = 10000;
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    int it = 0;
    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
	double c3 = c2*c2*x;
	c2 += d;
	c4 += d;
	a1 = c4 * a2 - c3 * a1;
	b1 = c4 * b2 - c3 * b1;

	c3 = c1 * c1 * x;
	c1 += d;
	c4 += d;
	a2 = c4 * a1 - c3 * a2;
	b2 = c4 * b1 - c3 * b2;

        if(trace) REprintf("it=%2d: ==> |b2|=%g", it, fabs(b2));
	if (fabs (b2) > scalefactor) {
            if(trace) REprintf("  Lrg |b2|");
	    a1 /= scalefactor;
	    b1 /= scalefactor;
	    a2 /= scalefactor;
	    b2 /= scalefactor;
	} else if (fabs (b2) < 1 / scalefactor) {
            if(trace) REprintf("  Sml |b2|");
	    a1 *= scalefactor;
	    b1 *= scalefactor;
	    a2 *= scalefactor;
	    b2 *= scalefactor;
	}
        if(trace) REprintf("\n");
        if(++it > maxit) {
            warning("non-convergence in logcf(), it = %d > maxit = %d", it, maxit);
	    break;
        }
    }
    if(trace && it <= maxit)
	REprintf("  logcf(*) used %d iterations.\n", it);
    return a2 / b2;
}


// To be  .Call()ed  from R :
SEXP R_logcf(SEXP x_, SEXP i_, SEXP d_, SEXP eps_, SEXP trace_)
{
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    // other args are checked via  asReal(.), asInteger() below
    R_xlen_t i,	n = XLENGTH(x_);
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    double *x = REAL(x_), *r = REAL(r_),
	ii = asReal(i_),
	d  = asReal(d_),
	eps= asReal(eps_);
    int trace = asInteger(trace_);
    if(ii <= 0.) error("i = %g <= 0", ii);
    if( d <  0.) error("d = %g <  0", d);

    for(i=0; i < n; i++) {
	r[i] = logcf(x[i], ii, d, eps, trace);
    }
    UNPROTECT(2);
    return r_;
}
