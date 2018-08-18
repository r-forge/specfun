#include "cdflib.h"

void V_cdfbet(int *which, double *p, double *q, double *x, double *y,
	      double *a, double *b,
	      int *status, double *bound, int *len);

void V_cdfbin(int *which, double *p, double *q, double *s, double *xn,
	      double *pr, double *ompr,
	      int *status, double *bound, int *len);

void V_cdfchi(int *which, double *p, double *q, double *x, double *df,
	      int *status, double *bound, int *len);

void V_cdfchn(int *which, double *p, double *q, double *x,
	      double *df, double *pnonc,
	      int *status, double *bound, int *len);

void V_cdff(int *which, double *p, double *q, double *f,
	    double *dfn, double *dfd,
	    int *status, double *bound, int *len);

void V_cdffnc(int *which, double *p, double *q, double *f,
	      double *dfn, double *dfd, double *phonc,
	      int *status, double *bound, int *len);

void V_cdfgam(int *which, double *p, double *q, double *x,
	      double *shape, double *scale,
	      int *status, double *bound, int *len);

void V_cdfnbn(int *which, double *p, double *q, double *s,
	      double *xn, double *pr, double *ompr,
	      int *status, double *bound, int *len);

void V_cdfnor(int *which, double *p, double *q, double *x,
	      double *mean, double *sd,
	      int *status, double *bound, int *len);

void V_cdfpoi(int *which, double *p, double *q, double *s, double *xlam,
	      int *status, double *bound, int *len);

void V_cdft(int *which, double *p, double *q, double *t, double *df,
	    int *status, double *bound, int *len);

void V_cdftnc(int *which, double *p, double *q, double *t,
	      double *df, double *pnonc,
	      int *status, double *bound, int *len);
