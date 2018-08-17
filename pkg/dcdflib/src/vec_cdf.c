/* ------------------------------------------------------------

 vectorized version of DCDFLIB  cdf**** functions.

 Use the same arguments as cdf***,  
 just an additional longint *len : simply loop through them *len times

 --- following the idea of  Donald H. MacQueen <macqueen1@llnl.gov>
 --- who did this for cdftnc only

 © 1999 Martin Maechler (C)
 This is GNU Public Licensed Software.
 ------------------------------------------------------------
*/

#include "vec_cdf.h"

void V_cdfbet(longint *which,double *p,double *q,double *x,double *y,
	      double *a,double *b,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfbet((int *)which,&p[i],&q[i],&x[i],&y[i],
	       &a[i],&b[i],(int *)&status[i],&bound[i]);
    }
}
 

void V_cdfbin(longint *which,double *p,double *q,double *s,double *xn,
	      double *pr,double *ompr,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfbin((int *)which,&p[i],&q[i],&s[i],&xn[i],
	       &pr[i],&ompr[i],(int *)&status[i],&bound[i]);
    }
}
 

void V_cdfchi(longint *which,double *p,double *q,double *x,double *df,
	      int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfchi((int *)which,&p[i],&q[i],&x[i],&df[i],
	       (int *)&status[i],&bound[i]);
    }
}
 
void V_cdfchn(longint *which,double *p,double *q,double *x,double *df,
	      double *pnonc,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfchn((int *)which,&p[i],&q[i],&x[i],&df[i],
	       &pnonc[i],(int *)&status[i],&bound[i]);
    }
}
 

void V_cdff(longint *which,double *p,double *q,double *f,double *dfn,
	    double *dfd,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdff((int *)which,&p[i],&q[i],&f[i],&dfn[i],
	     &dfd[i],(int *)&status[i],&bound[i]);
    }
}
 
void V_cdffnc(longint *which,double *p,double *q,double *f,double *dfn,
	      double *dfd,double *phonc,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdffnc((int *)which,&p[i],&q[i],&f[i],&dfn[i],
	       &dfd[i],&phonc[i],(int *)&status[i],&bound[i]);
    }
}
 

void V_cdfgam(longint *which,double *p,double *q,double *x,double *shape,
	      double *scale,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfgam((int *)which,&p[i],&q[i],&x[i],&shape[i],
	       &scale[i],(int *)&status[i],&bound[i]);
    }
}
 

void V_cdfnbn(longint *which,double *p,double *q,double *s,double *xn,
	      double *pr,double *ompr,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfnbn((int *)which,&p[i],&q[i],&s[i],&xn[i],
	       &pr[i],&ompr[i],(int *)&status[i],&bound[i]);
    }
}
 

void V_cdfnor(longint *which,double *p,double *q,double *x,double *mean,
	      double *sd,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfnor((int *)which,&p[i],&q[i],&x[i],&mean[i],
	       &sd[i],(int *)&status[i],&bound[i]);
    }
}
 

void V_cdfpoi(longint *which,double *p,double *q,double *s,double *xlam,
	      int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdfpoi((int *)which,&p[i],&q[i],&s[i],&xlam[i],
	       (int *)&status[i],&bound[i]);
    }
}
 

void V_cdft(longint *which,double *p,double *q,double *t,double *df,
	    int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdft((int *)which,&p[i],&q[i],&t[i],&df[i],
	     (int *)&status[i],&bound[i]);
    }
}
 
void V_cdftnc(longint *which,double *p,double *q,double *t,double *df,
	      double *pnonc,int *status,double *bound, longint *len)
{
    int i;
    for (i = 0; i < *len; i++) {
	cdftnc((int *)which,&p[i],&q[i],&t[i],&df[i],
	       &pnonc[i],(int *)&status[i],&bound[i]);
    }
}
