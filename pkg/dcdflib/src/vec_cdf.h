#include "cdflib.h"
/* #include "cdflib-new.h"
 *                 ^^^^  only momentarily:*/

void V_cdfbet(longint *which,double *p,double *q,double *x,double *y,
	      double *a,double *b,
	      int *status,double *bound, longint *len);

void V_cdfbin(longint *which,double *p,double *q,double *s,double *xn,
	      double *pr,double *ompr,
	      int *status,double *bound, longint *len);

void V_cdfchi(longint *which,double *p,double *q,double *x,double *df,
	      int *status,double *bound, longint *len);

void V_cdfchn(longint *which,double *p,double *q,double *x,
	      double *df, double *pnonc,
	      int *status,double *bound, longint *len);

void V_cdff(longint *which,double *p,double *q,double *f,
	    double *dfn, double *dfd,
	    int *status,double *bound, longint *len);

void V_cdffnc(longint *which,double *p,double *q,double *f,
	      double *dfn, double *dfd,double *phonc,
	      int *status,double *bound, longint *len);

void V_cdfgam(longint *which,double *p,double *q,double *x,
	      double *shape, double *scale,
	      int *status,double *bound, longint *len);

void V_cdfnbn(longint *which,double *p,double *q,double *s,
	      double *xn, double *pr,double *ompr,
	      int *status,double *bound, longint *len);

void V_cdfnor(longint *which,double *p,double *q,double *x,
	      double *mean, double *sd,
	      int *status,double *bound, longint *len);

void V_cdfpoi(longint *which,double *p,double *q,double *s, double *xlam,
	      int *status,double *bound, longint *len);

void V_cdft(longint *which,double *p,double *q,double *t,double *df,
	    int *status,double *bound, longint *len);

void V_cdftnc(longint *which,double *p,double *q,double *t,
	      double *df, double *pnonc,
	      int *status,double *bound, longint *len);
