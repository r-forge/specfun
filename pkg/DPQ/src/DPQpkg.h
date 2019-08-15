
#include <Rmath.h>
// for F77_NAME() :
#include <R_ext/RS.h>

// for SEXP:
#include <Rinternals.h>

// From R's source
// -------- excerpt from  nmath.h ---------------------------------
#include <R_ext/Error.h>
# define MATHLIB_ERROR(fmt,x)		error(fmt,x);
# define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) warning(fmt,x,x2,x3,x4,x5)

#include <R_ext/Arith.h>
#define ML_POSINF	R_PosInf
#define ML_NEGINF	R_NegInf
#define ML_NAN		R_NaN

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) gettext (String)
#else
#define _(String) (String)
#endif
// ----------- end of excerpt from R's  nmath.h -------------------------

#include "dpq.h"
//        =====

/* Type 'logical':
   1) use for .C() called functions for clarity
   2) remember which were Fortran logicals (now in C)
*/
typedef int logical;


// qchisq_appr.c : -------------------------------------------------------------

void qchisq_appr_v(double *P, int *n, double *nu, double *tol,
		   logical *lower_tail, logical *log_p,
		   /* result: */ double *q)
    ;

// pnchisq-it.c : --------------------------------------------------------------

void Pnchisq_it(double *x, double *f, double *theta,
		/* FIXME?? additionally vectorize in (x, prob) or (x,f,th) |-> prob */
		double *errmax, double *reltol, int *itrmax,
		int *i_0, int *n_terms, double *terms, double *prob)
    ;

// 310-pnbeta.c : --------------------------------------------------------------

void ncbeta(double *a, double *b, double *lambda, double *x, int *n,
	    double *errmax, int *itrmax, int *ifault, double *res)
    ;

// ppois-direct.c : ------------------------------------------------------------

void ppois_D(double *X, int *n, double *lambda, /* result: */ double *prob)
    ;

// wienergerm_nchisq.c : -------------------------------------------------------

// TODO: export h() function ? (with longer name)?

double nonc_chi(double x, double ncp, double df, int lower_tail, int log_p,
		int variant);
// Called via .C():
void pchisqV(double *x, int *n, /* vectorized in x : x[1..n] : */
	     double *ncp, double *df,
	     logical *lower_tail, logical *log_p, int *variant)
    ;

// wienergerm_nchisq_F.f : -----------------------------------------------------
int F77_NAME(noncechi)(int *variant,
		       double *argument, double *noncentr, double *df, double *p,
		       int *ifault);


// algdiv.c: --------------------------------------------------------------------

double algdiv(double a, double b);
// .Call()ed :
SEXP R_algdiv(SEXP a_, SEXP b_)
    ;
