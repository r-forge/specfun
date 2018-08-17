/* -- $Id: cdflib.h,v 1.12 2018/08/17 13:42:43 maechler Exp $ */

/* TODO:

   o replaced 'unsigned long' by 'long'

   o Add 'void' [R/S interface] functions also for things like  algdiv()

   o add comments here


   o find calling structure
 */

/* NB:  V_cdftnc(..) etc are all in
 * --  ./vec_cdf.h
 *       ~~~~~~~~~ */


/* The following is from R's  "S.h" -- S/R compatibility : */
#ifndef longint
#define longint int
#endif

double algdiv(double*,double*);
double alngam(double*);
double alnrel(double*);
double apser(double*,double*,double*,double*);
double basym(double*,double*,double*,double*);
double bcorr(double*,double*);
double betaln(double*,double*);
double bfrac(double*,double*,double*,double*,double*,double*);
double bpser(double*,double*,double*,double*);
double brcmp1(int*,double*,double*,double*,double*);
double brcomp(double*,double*,double*,double*);
double bup(double*,double*,double*,double*,int*,double*);

void bgrat (double*,double*,double*,double*,double*,double*,int*);
void bratio(double*,double*,double*,double*,double*,double*,int*);

void cdfbet(int*,double*,double*,double*,double*,double*,double*,int*,double*);
void cdfbin(int*,double*,double*,double*,double*,double*,double*,int*,double*);
void cdfchi(int*,double*,double*,double*,double*,        int*,double*);
void cdfchn(int*,double*,double*,double*,double*,double*,int*,double*);
void cdff  (int*,double*,double*,double*,double*,double*,        int*,double*);
void cdffnc(int*,double*,double*,double*,double*,double*,double*,int*,double*);
void cdfgam(int*,double*,double*,double*,double*,double*,int*,double*);
void cdfnbn(int*,double*,double*,double*,double*,double*,double*,int*,double*);
void cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
void cdfpoi(int*,double*,double*,double*,double*,int*,double*);
void cdft  (int*,double*,double*,double*,double*,        int*,double*);
void cdftnc(int*,double*,double*,double*,double*,double*,int*,double*);

void cumbet(double*,double*,double*,double*,double*,double*);
void cumbin(double*,double*,double*,double*,double*,double*);
void cumchi(double*,double*,double*,double*);
void cumchn(double*,double*,double*,double*,double*);
void cumf  (double*,double*,double*,double*,double*);
void cumfnc(double*,double*,double*,double*,double*,double*);
void cumgam(double*,double*,double*,double*);
void cumnbn(double*,double*,double*,double*,double*,double*);
void cumnor(double*,double*,double*);
void cumpoi(double*,double*,double*,double*);
void cumt  (double*,double*,double*,double*);
void cumtnc(double*,double*,double*,double*,double*);

double devlpl(double*, int,double);
double dinvnr(double *p,double *q);
void dinvr   (int*,double*,double*,unsigned long*,unsigned long*);
void dstinv  (double*,double*,double*,double*,double*,double*,double*);
double dt1(double*,double*,double*);
/*
static void E0000(int,int*,double*,double*,unsigned long*,
	   unsigned long*,double*,double*,double*,
	   double*,double*,double*,double*);
static void E0001(int,int*,double*,double*,double*,double*,
	   unsigned long*,unsigned long*,double*,double*,
	   double*,double*);
*/
void dzror(int*,double*,double*,double*,double *,
           unsigned long*,unsigned long*);
void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl);

double erf1(double*);
double erfc1(int*,double*);
double erfi(double*, double*);
double esum(int*,double*);
double fpser(double*,double*,double*,double*);
double gam1(double*);
void gaminv(double*,double*,double*,double*,double*,int*);
double gamln(double*);
double gamln1(double*);
double Xgamm(double*);
void grat1(double*,double*,double*,double*,double*,double*);
void gratio(double*,double*,double*,double*,int*);
double gsumln(double*,double*);
void pni(double *, double *, double *, double *, long *);
double psi(double*);
double rcomp(double*,double*);
double cdf_rexp(double*);
double rlog(double*);
double rlog1(double*);
double stvaln(double*);

double fifdint(double);
double fifdmax1(double,double);
double fifdmin1(double,double);
double fifdsign(double,double);
long fifidint(double);
long fifmod(long,long);
#ifdef _USE_FTNSTOP_
void ftnstop(char*);
#endif

double exparg(int);
double spmpar(int);
extern int ipmpar(int);

/*----------------- NEW (M.Maechler) ----------- */

#ifdef DEBUG
# define DBGprt1(ch) printf(ch)
# define DBGprt2(ch,x) printf(ch,x)
# define DBGprt3(ch,x,y) printf(ch,x,y)
# define DBGprt4(ch,x,y,z) printf(ch,x,y,z)
# define DBGprt5(ch,w,x,y,z) printf(ch,w,x,y,z)
# define DBGprt6(ch,v,w,x,y,z) printf(ch,v,w,x,y,z)
#else
# define DBGprt1(ch)
# define DBGprt2(ch,x)
# define DBGprt3(ch,x,y)
# define DBGprt4(ch,x,y,z)
# define DBGprt5(ch,w,x,y,z)
# define DBGprt6(ch,v,w,x,y,z)
#endif

/* From   Mathlib.h */
#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif
#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

