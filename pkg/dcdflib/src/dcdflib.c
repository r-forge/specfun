/* -- $Id: dcdflib.c,v 1.35 2018/08/17 13:43:48 maechler Exp $
 *
 * The original was in Barry Brown's
 *   ftp://odin.mdacc.tmc.edu/pub/source/dcdflib.c-1.1.tar.gz
 *   (dated Dec 8, 1997; length 61435)
 *
 * >>> Call Tree -- via 'ftncheck' from original Fortran sources
 * >>> /u/sfs/div/source/BBrown/dcdflib.f/src/All-calltree
 *
 * Martin Maechler: I found that the UNDERLYING NSWC Library
 * ---------------  has had revisions that were not yet incorporated in
 * dcdlib 1.1
 * ---------- > /u/sfs/div/source/BBrown/dcdflib.f/src/diff-NSWC.out   and
 *		/usr/local/lib.src/NSWC/nswc/
 */

/* remember which were Fortran logicals : */
#define logical unsigned long

#include <R.h>
#include <Rmath.h>

#include "cdflib.h"
/* now includes the  DBGprt<n> macros */

/*
 * A comment about ints and longs - whether ints or longs are used should
 * make no difference, but where double r-values are assigned to ints the
 * r-value is cast converted to a long, which is then assigned to the int
 * to be compatible with the operation of fifidint.
 */

/* Global Constants */

/* LG_c. : Minimax coefficients used for DEL(x) in  ln(gamma(x)) for larger x */
static double LG_c0 =  .833333333333333e-01;
static double LG_c1 = -.277777777760991e-02;
static double LG_c2 =  .793650666825390e-03;
static double LG_c3 = -.595202931351870e-03;
static double LG_c4 =  .837308034031215e-03;
static double LG_c5 = -.165322962780713e-02;

/* when pointers must be passed... */
static int I0 = 0;
static int I1 = 1;


double algdiv(double *a,double *b)
/*
-----------------------------------------------------------------------

     Computation of ln(gamma(b)/gamma(a+b)) when   b >= 8

			 --------

     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
	{R's Mathlib calls DEL(X) === lgammacor(x)}
-----------------------------------------------------------------------
*/
{
  /* NB: This is VERY similar to  bcorr() below !!! ----*/

    static double c,d,h,s11,s3,s5,s7,s9,t,u,v,w,x,x2,T1;
/*
     ..
     .. Executable Statements ..
*/
    if(*a > *b) {
	h = *b/ *a;
	c = 1/(1+h);
	x = h/(1+h);
	d = *a+(*b-0.5);
    }
    else {
	h = *a/ *b;
	c = h/(1+h);
	x = 1/(1+h);
	d = *b+(*a-0.5);
    }

/*		  SET SN = (1 - X^N)/(1 - X)  */
    x2 = x*x;
    s3 = 1+(x+x2);
    s5 = 1+(x+x2*s3);
    s7 = 1+(x+x2*s5);
    s9 = 1+(x+x2*s7);
    s11= 1+(x+x2*s9);
/*
		SET W = DEL(B) - DEL(A + B)
*/
    t = pow(1/ *b,2);
    w = ((((LG_c5*s11*t+LG_c4*s9)*t+LG_c3*s7)*t+LG_c2*s5)*t+LG_c1*s3)*t+LG_c0;
    w *= (c/ *b);
/*
		    COMBINE THE RESULTS
*/
    T1 = *a/ *b;
    u = d*alnrel(&T1);
    v = *a*(log(*b)-1);
    if(u > v)
	return w-v-u;
    else
	return w-u-v;
}

double bcorr(double *a0,double *b0)
{
/*
-----------------------------------------------------------------------

     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8.
	{R's Mathlib calls DEL(X) === lgammacor(x)}

-----------------------------------------------------------------------
 NB:  Similar to  algdiv()  above !
*/
static double a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;

    a = fmin2(*a0,*b0);
    b = fmax2(*a0,*b0);
    h = a/b;
    c = h/(1+h);
    x = 1/(1+h);
/*
		SET SN = (1 - X^N)/(1 - X)
*/
    x2 = x*x;
    s3 = 1+(x+x2);
    s5 = 1+(x+x2*s3);
    s7 = 1+(x+x2*s5);
    s9 = 1+(x+x2*s7);
    s11= 1+(x+x2*s9);
/*
		SET W = DEL(B) - DEL(A + B)
*/
    t = pow(1/b,2);
    w = ((((LG_c5*s11*t+LG_c4*s9)*t+LG_c3*s7)*t+LG_c2*s5)*t+LG_c1*s3)*t+LG_c0;
    w *= (c/b);
/*
		COMPUTE	 DEL(A) + W
*/
    t = pow(1/a,2);
    return (((((LG_c5*t+LG_c4)*t+LG_c3)*t+LG_c2)*t+LG_c1)*t+LG_c0)/a + w;
} /* bcorr() */


double alngam(double *x)
{
/*
**********************************************************************

     double precision LN of the GAMma function


			      Function


     Returns the natural logarithm of GAMMA(X).


			      Arguments


     X --> value at which scaled log gamma is to be returned
		    X is DOUBLE PRECISION


			      Method


     If X <= 6, then use recursion to get X below 3
     then apply rational approximation number 5236 of
     Hart et al, Computer Approximations, John Wiley and
     Sons, NY, 1968.

     If X > 6, then use recursion to get X to at least 12 and
     then use formula 5423 of the same source.

**********************************************************************
*/
#define hln2pi 0.91893853320467274178
    static double coef[5] = {
	0.83333333333333023564e-1,-0.27777777768818808e-2,0.79365006754279e-3,
	-0.594997310889e-3,0.8065880899e-3
    };
    static double scoefd[4] = {
	62.003838007126989331,9.822521104713994894,-8.906016659497461257, 1.
    };
    static double scoefn[9] = {
	62.003838007127258804,36.036772530024836321,20.782472531792126786,
	 6.338067999387272343, 2.15994312846059073, 0.3980671310203570498,
	0.1093115956710439502,.0092381945590275995, .0029737866448101651
    };
    double alngam,offset,prod,xx;
    int i,n;
/*
  ..
  .. Executable Statements ..
*/
    if(*x <= 6) {
	prod = 1;
	xx = *x;
	while(xx > 3) {
	    xx -= 1;
	    prod *= xx;
	}

	if(*x < 2)
	    while(xx < 2) {
		prod /= xx;
		xx += 1;
	    }

	alngam = devlpl(scoefn,9,xx-2) / devlpl(scoefd,4,xx-2);
	/*
	  COMPUTE RATIONAL APPROXIMATION TO GAMMA(X)
	*/
	alngam *= prod;
	alngam = log(alngam);
    }
    else { /*  x > 6 */

	offset = hln2pi;
	/*
	  IF NECESSARY MAKE X AT LEAST 12 AND CARRY CORRECTION IN OFFSET
	*/
	n = (long)(12-*x);
	if(n > 0) {
	    for(i=0, prod=1; i<n; i++)
		prod *= (*x+ i);
	    offset -= log(prod);
	    xx = *x + n;
	} else {
	    xx = *x;
	}
	/*
	  COMPUTE POWER SERIES
	*/
	alngam = devlpl(coef, 5, 1/pow(xx,2))/xx;
	alngam += (offset+(xx-0.5)*log(xx)-xx);
    }
    return alngam;
#undef hln2pi
}

double alnrel(double *a)
/*
-----------------------------------------------------------------------
	    LN(1 + A)
-----------------------------------------------------------------------
*/
{
static double p1 = -1.29418923021993;
static double p2 =  .405303492862024;
static double p3 = -.0178874546012214;
static double q1 = -1.62752256355323;
static double q2 =  .747811014037616;
static double q3 = -.0845104217945565;

    double t,t2,w;

    if(fabs(*a) <= 0.375) {
	t = *a/(*a+2);
	t2 = t*t;
	w = (((p3*t2+p2)*t2+p1)*t2+1)/
	    (((q3*t2+q2)*t2+q1)*t2+1);
	return 2*t*w;
    } else {
	if(*a < 0) t = (*a + .5) + .5; else t = 1. + *a;
	return log(t);
    }
}

double apser(double *a,double *b,double *x,double *eps)
/*
-----------------------------------------------------------------------
     INCOMPLETE BETA RATIO  I_{1-X}(B,A)  when A is VERY SMALL.

     Use only if
		A <= min(EPS,EPS*B),	B*X <= 1,  and	X <= 0.5.
-----------------------------------------------------------------------
*/
{
    static double g = .577215664901533;
    double aj,bx,c,s,t,tol;
    int j;

    bx = *b * *x;
    t = *x-bx;
    if(*b * *eps <= 0.02)
	c = log(*x)+psi(b)+g+t;
    else
	c = log(bx)+g+t;

    tol = 5**eps*fabs(c);
    j = 1;
    s = 0;
    do {
	j += 1;
	t *= (*x-bx/j);
	aj = t/j;
	s += aj;
    } while(fabs(aj) > tol);

    return -(*a*(c+s));
}

double basym(double *a,double *b,double *lambda,double *eps)
/*
-----------------------------------------------------------------------
     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.

     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
     A AND B ARE GREATER THAN OR EQUAL TO 15.
-----------------------------------------------------------------------
*/
{
    /* Constants */
    static double e0 = 1.12837916709551;/* = 2/SQRT(PI) */
    static double e1 = .353553390593274;/* = 2^(-3/2)	*/
    static int one = 1;

#define num 20/* num := maximum value that  n  can take in the outermost
  for() loop below. It is required that num be even.
  The arrays a0, b0, c, d have dimension num + 1.
*/
    double a0[21],b0[21],c[21],d[21];

    double bsum,dsum,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,t0,t1,w,w0,z,z0,z2,zn,znm1,
	T1,T2;
    int i,imj,j,m,mmj,n;

    if(*a < *b) {
	h = *a/ *b;
	r0 = 1/(1+h);
	r1 = (*b-*a)/ *b;
	w0 = 1/sqrt(*a*(1+h));
    } else {
	h = *b/ *a;
	r0 = 1/(1+h);
	r1 = (*b-*a)/ *a;
	w0 = 1/sqrt(*b*(1+h));
    }

    T1 = -(*lambda/ *a);
    T2 = *lambda/ *b;
    f = *a*rlog1(&T1)+*b*rlog1(&T2);
    t = exp(-f);
    if(t == 0) return t;
    z0 = sqrt(f);
    z = 0.5*(z0/e1);
    z2 = f+f;
    a0[0] = 2/3*r1;
    c[0] = -(0.5*a0[0]);
    d[0] = -c[0];
    j0 = 0.5/e0*erfc1(&one,&z0);
    j1 = e1;
    sum = j0+d[0]*w0*j1;
    s = 1;
    h2 = h*h;
    hn = 1;
    w = w0;
    znm1 = z;
    zn = z2;
    for(n=2; n<=num; n+=2) {
	hn *= h2;
	a0[n-1] = 2*r0*(1+h*hn)/(n+2);
	s += hn;
	a0[n] = 2*r1*s/(n+3);
	for(i=n; i<=n+1; i++) { /* two times */
	    r = -(0.5*(i+1));
	    b0[0] = r*a0[0];
	    for(m=2; m<=i; m++) {
		bsum = 0;
		for(j=1; j<m; j++) {
		    mmj = m-j;
		    bsum += ((j*r-mmj)*a0[j-1]*b0[mmj-1]);
		}
		b0[m-1] = r*a0[m-1] + bsum/m;
	    }
	    c[i-1] = b0[i-1]/(i+1);
	    dsum = 0;
	    for(j=1; j<i; j++) {
		imj = i-j;
		dsum += (d[imj-1]*c[j-1]);
	    }
	    d[i-1] = -(dsum+c[i-1]);
	}
	j0 = e1*znm1+ ((double)n-1)*j0;
	j1 = e1*zn  + ((double) n )*j1;
	znm1 *= z2;
	zn *= z2;
	w *= w0;
	t0 = d[n-1]*w*j0;
	w *= w0;
	t1 = d[n]*w*j1;
	sum += (t0+t1);
	if(fabs(t0)+fabs(t1) <= *eps*sum) break;
    }

    return e0*t* exp(-bcorr(a,b)) * sum;
} /* basym() */


double betaln(double *a0,double *b0)
/*
-----------------------------------------------------------------------
     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
-----------------------------------------------------------------------
*/
{
static double e = .918938533204673;/* E := 0.5*LN(2*PI) */

    double a,b,c,h,u,v,w,z;
    int i,n;
    double T1;

    a = fmin2(*a0,*b0);
    b = fmax2(*a0,*b0);

    if(a < 8) {
	if(a < 1) {
	    /*-----------------------------------------------------------
	      PROCEDURE WHEN A < 1
	      ----------------------------------------------------------- */
	    if(b < 8) {
		T1 = a+b;
		return gamln(&a)+(gamln(&b)-gamln(&T1));
	    }
	    return gamln(&a)+algdiv(&a,&b);
	}
	/*---------------------------------------------------------------
	  PROCEDURE WHEN 1 <= A < 8
	  --------------------------------------------------------------- */
	if(a <= 2) {
	    if(b <= 2)
		return gamln(&a)+gamln(&b)-gsumln(&a,&b);
	    /* b > 2: */
		w = 0;
		if(b < 8) goto S60;
		return gamln(&a)+algdiv(&a,&b);
	}

	/*
	  REDUCTION OF A in (2,8)  WHEN B <= 1000
	*/
	if(b <= 1000) {
	    n = (long)(a - 1);
	    w = 1;
	    for(i=1; i<=n; i++) {
		a -= 1;
		h = a/b;
	    w *= (h/(1+h));
	    }
	    w = log(w);
	    if(b >= 8)
		return w+gamln(&a)+algdiv(&a,&b);


	S60:
	    /*
	      REDUCTION OF B WHEN B < 8
	    */
	    n = (long)(b - 1);
	    z = 1;
	    for(i=1; i<=n; i++) {
		b -= 1;
		z *= (b/(a+b));
	    }
	    return w+log(z)+(gamln(&a)+(gamln(&b)-gsumln(&a,&b)));
	}

	/*
	  REDUCTION OF A WHEN B > 1000
	*/
	n = (long)(a - 1);
	w = 1;
	for(i=1; i<=n; i++) {
	    a -= 1;
	    w *= (a/(1+a/b));
	}
	return	log(w)-(double)n*log(b)+(gamln(&a)+algdiv(&a,&b));
    }

    else {
/*
  -----------------------------------------------------------------------
  PROCEDURE WHEN A >= 8
  -----------------------------------------------------------------------*/
	w = bcorr(&a,&b);
	h = a/b;
	c = h/(1+h);
	u = -((a-0.5)*log(c));
	v = b*alnrel(&h);
	if(u > v)
	    return -(0.5*log(b))+e+w-v-u;
	else
	    return -(0.5*log(b))+e+w-u-v;
    }
} /* betaln() */


double bfrac(double *a,double *b,double *x,double *y,double *lambda,
	     double *eps)
{
/*
-----------------------------------------------------------------------
     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B > 1.
     IT IS ASSUMED THAT	 LAMBDA = (A + B)*Y - B.
-----------------------------------------------------------------------
*/
    double bfrac,alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,t,w,yp1;

    bfrac = brcomp(a,b,x,y);
    if(bfrac == 0) return bfrac;
    c = 1+*lambda;
    c0 = *b/ *a;
    c1 = 1+1/ *a;
    yp1 = *y+1;
    n = 0;
    p = 1;
    s = *a+1;
    an = 0;
    bn = anp1 = 1;
    bnp1 = c/c1;
    r = c1/c;
/*
	CONTINUED FRACTION CALCULATION
*/
    do {
	n += 1;
	t = n/ *a;
	w = n*(*b-n)**x;
	e = *a/s;
	alpha = p*(p+c0)*e*e*(w**x);
	if(alpha <= 0) break;/* << NEW from NSWC; MM */
	e = (1+t)/(c1+t+t);
	beta = n+w/s+e*(c+n*yp1);
	p = 1+t;
	s += 2;
	/*
	  UPDATE AN, BN, ANP1, AND BNP1
	*/
	t = alpha*an+beta*anp1;
	an = anp1;
	anp1 = t;
	t = alpha*bn+beta*bnp1;
	bn = bnp1;
	bnp1 = t;
	r0 = r;
	r = anp1/bnp1;
	if(fabs(r-r0) <= *eps*r) break;
	/*
	  Rescale  an, bn, anp1 & bnp1
	*/
	an /= bnp1;
	bn /= bnp1;
	anp1 = r;
	bnp1 = 1;
    } while(1);
/*
		 TERMINATION
*/
    bfrac *= r;
    return bfrac;
} /* bfrac() */

void bgrat(double *a,double *b,double *x,double *y,double *w,
	   double *eps,int *ierr)
{
/*
-----------------------------------------------------------------------
     Asymptotic Expansion for IX(A,B) when  A >= 15 > 1 >= B.

     The result of the expansion is added to W.
     EPS is the tolerance used.
     IERR is a variable that reports the status of the results.
-----------------------------------------------------------------------
*/
    double c[30],d[30];

    double bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,t2,u,v,z, T1, tol;
    int i,n;

    bm1 = *b-0.5-0.5;
    nu = *a+0.5*bm1;
    if(*y > 0.375) {
	lnx = log(*x);
    } else {
	T1 = -*y;
	lnx = alnrel(&T1);
    }

    z = -(nu*lnx);
    if(*b*z == 0) goto L_error;
/*
		COMPUTATION OF THE EXPANSION
		SET R = EXP(-Z)* Z^B / GAMMA(B)
*/
    r = *b*(1+gam1(b))*exp(*b*log(z));
    r *= (exp(*a*lnx)*exp(0.5*bm1*lnx));
    u = algdiv(b,a) + *b*log(nu);
    u = r*exp(-u);
    if(u == 0) goto L_error;
    grat1(b,&z,&r,&p,&q,eps);

    tol = 15* *eps; /* << NEW from NSWC	 (and use 'tol', not *eps below; MM */
    v = 0.25 /(nu*nu);
    t2 = 0.25*lnx*lnx;
    l = *w/u;
    j = q/r;
    sum = j;
    t = cn = 1;
    n2 = 0;
    for(n=1; n <= 30; n++) {
	bp2n = *b+n2;
	j = (bp2n*(bp2n+1)*j+(z+bp2n+1)*t)*v;
	n2 += 2;
	t *= t2; /* = t2 ^ n =	(lnx/2) ^ (2n) */
	cn /= (n2*(n2+1));
	c[n-1] = cn;
	s = 0;
	if(n != 1) {
	    coef = *b-(double)n;
	    for(i=1; i < n; i++) {
		s += (coef*c[i-1]*d[n-i-1]);
		coef += *b;
	    }
	}
	d[n-1] = bm1*cn+s/(double)n;
	dj = d[n-1]*j;
	sum += dj;
	if(sum <= 0) goto L_error;
	if(fabs(dj) <= tol*(sum+l)) break;
    }
/*
		ADD THE RESULTS TO W
*/
    *ierr = 0;
    *w += (u*sum);
    return;

L_error:/*	THE EXPANSION CANNOT BE COMPUTED */
    *ierr = 1;
    return;

} /* bgrat() */

double bpser(double *a,double *b,double *x,double *eps)
{
/*
-----------------------------------------------------------------------
     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B <= 1
     OR B*X <= 0.7.  EPS IS THE TOLERANCE USED.
-----------------------------------------------------------------------
*/
    double bpser,a0,apb,b0,c,n,sum,t,tol,u,w,z;
    int i,m;

    if(*x == 0) return 0;

    bpser = 0;
 /*-----------------------------------------------------------------------
	    COMPUTE THE FACTOR X^A/(A*BETA(A,B))
  -----------------------------------------------------------------------*/
    a0 = fmin2(*a,*b);
    if(a0 >= 1) {
	z = *a*log(*x)-betaln(a,b);
	bpser = exp(z)/ *a;
    }
    else {

	b0 = fmax2(*a,*b);

	if(b0 < 8) {
	    if(b0 <= 1) {/*	procedure for A0 < 1 and B0 <= 1 */

		bpser = pow(*x,*a);
		if(bpser == 0) return bpser;
		apb = *a+*b;
		if(apb <= 1) {
		    z = 1+gam1(&apb);
		}
		else {
		    u = *a+*b-1.e0;
		    z = (1+gam1(&u))/apb;
		}

		c = (1+gam1(a))*(1+gam1(b))/z;
		bpser *= (c*(*b/apb));
	    }
	    else {/*	procedure for A0 < 1 and 1 < B0 < 8 */

		u = gamln1(&a0);
		m = (long)(b0 - 1);
		if(m >= 1) {
		    c = 1;
		    for(i=1; i<=m; i++) {
			b0 -= 1;
			c *= (b0/(a0+b0));
		    }
		    u = log(c)+u;
		}
		z = *a*log(*x)-u;
		b0 -= 1;
		apb = a0+b0;
		if(apb <= 1) {
		    t = 1+gam1(&apb);
		}
		else {
		    u = a0+b0-1.e0;
		    t = (1+gam1(&u))/apb;
		}

		bpser = exp(z)*(a0/ *a)*(1+gam1(&b0))/t;
	    }
	}
	else {/*	procedure for A0 < 1 and B0 >= 8 */

	    u = gamln1(&a0)+algdiv(&a0,&b0);
	    z = *a*log(*x)-u;
	    bpser = a0/ *a*exp(z);
	}
    }

    if(bpser == 0 || *a <= 0.1**eps) return bpser;
/*
-----------------------------------------------------------------------
		     COMPUTE THE SERIES
-----------------------------------------------------------------------
*/
    sum = n = 0;
    c = 1;
    tol = *eps/ *a;
    do {
	n += 1;
	c *= ((0.5+(0.5-*b/n))**x);
	w = c/(*a+n);
	sum += w;
    } while(fabs(w) > tol);
    bpser *= (1+*a*sum);
    return bpser;
} /* bpser() */

void bratio(double *a,double *b, double *x,double *y,
	    double *w,double *w1, int *ierr)
{
/* S or R :  w = pbeta(x, a,b)
-----------------------------------------------------------------------

	    EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)

		     --------------------

     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X <= 1
     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES

		      W	 = IX(A,B)
		      W1 = 1 - IX(A,B)

     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
     ONE OF THE FOLLOWING VALUES ...

	IERR = 1  IF A OR B IS NEGATIVE
	IERR = 2  IF A = B = 0
	IERR = 3  IF X < 0 OR X > 1
	IERR = 4  IF Y < 0 OR Y > 1
	IERR = 5  IF X + Y != 1
	IERR = 6  IF X = A = 0
	IERR = 7  IF Y = B = 0

--------------------
     WRITTEN BY ALFRED H. MORRIS, JR.
	NAVAL SURFACE WARFARE CENTER
	DAHLGREN, VIRGINIA
     REVISED ... NOV 1991
-----------------------------------------------------------------------
*/
    int ierr1,ind,n;
    double a0,b0,eps,lambda,t,x0,y0,z;
    double T2,T3,T4,T5;

/*
     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
	    FLOATING POINT NUMBER FOR WHICH 1 + EPS > 1
*/
    eps = spmpar(1);
    /*
    DBGprt4("N> bratio(%g,%g, %g) ", *a, *b, *x);
    */
    *w = *w1 = 0;

    /* Argument Checking & Error Returns : */
    if(*a < 0 || *b < 0) {	*ierr = 1; return;}
    if(*a == 0 && *b == 0) {	*ierr = 2; return;}
    if(*x < 0 || *x > 1) {	*ierr = 3; return;}
    if(*y < 0 || *y > 1) {	*ierr = 4; return;}
    z = *x+*y-0.5-0.5;
    if(fabs(z) > 3*eps) {	*ierr = 5; return;}

    *ierr = 0;
    if(*x == 0) goto x_is_0;
    if(*y == 0) goto y_is_0;
    if(*a == 0) goto a_is_0;
    if(*b == 0) goto b_is_0;
    eps = fmax2(eps,1.e-15);
    if(fmax2(*a,*b) < 1.e-3*eps) goto L_ab_small;
    ind = 0;
    a0 = *a;
    b0 = *b;
    x0 = *x;
    y0 = *y;
    if(a0 <= 1 || b0 <= 1) {/* procedure for  A0 <= 1 or  B0 <= 1 */
	if(*x > 0.5) {
	    ind = 1;
	    a0 = *b;
	    b0 = *a;
	    x0 = *y;
	    y0 = *x;
	}
	if(b0 < fmin2(eps,eps*a0)) {		/* S90 */
	    *w = fpser(&a0,&b0,&x0,&eps);
	    *w1 = 0.5+(0.5-*w);
	    goto End_ok;
	}
	if(a0 < fmin2(eps,eps*b0) && b0*x0 <= 1) { /* S100 */
	    *w1 = apser(&a0,&b0,&x0,&eps);
	    *w = 0.5+(0.5-*w1);
	    goto End_ok;
	}
	if(fmax2(a0,b0) <= 1) {
	    if(a0 >= fmin2(0.2,b0)) goto S110;
	    if(pow(x0,a0) <= 0.9) goto S110;
	    if(x0 >= 0.3) goto S120;
	}
	else {
	    if(b0 <= 1) goto S110;
	    if(x0 >= 0.3) goto S120;
	    if(x0 < 0.1) {
		if(pow(x0*b0,a0) <= 0.7) goto S110;
	    }

	    if(b0 > 15) goto S150;
	}
	n = 20;
	goto S140;
    }
    else {/*		procedure for A0 > 1 AND B0 > 1 */
	if(*a <= *b)
	    lambda = *a-(*a+*b)**x;
	else
	    lambda = (*a+*b)**y-*b;

	if(lambda < 0) {
	    ind = 1;
	    a0 = *b;
	    b0 = *a;
	    x0 = *y;
	    y0 = *x;
	    lambda = fabs(lambda);
	}

	if(b0 < 40 && b0*x0 <= 0.7) goto S110;
	if(b0 < 40) goto S160;
	if(a0 > b0) { /* S80 */
	    if(b0 <= 100 || lambda > 0.03*b0) goto S130;
	    goto Use_basym;
	}
	else {
	    if(a0 <= 100 || lambda > 0.03*a0) goto S130;
	    goto Use_basym;
	}
    }

    /*
      EVALUATION OF THE APPROPRIATE ALGORITHM
    */
  S110:
    *w = bpser(&a0,&b0,&x0,&eps);
    *w1 = 0.5+(0.5-*w);
    goto End_ok;
  S120:
    *w1 = bpser(&b0,&a0,&y0,&eps);
    *w = 0.5+(0.5-*w1);
    goto End_ok;
  S130:
    T2 = 15*eps;
    *w = bfrac(&a0,&b0,&x0,&y0,&lambda,&T2);
    *w1 = 0.5+(0.5-*w);
    goto End_ok;

  S140:
    *w1 = bup(&b0,&a0,&y0,&x0,&n,&eps);
    b0 += (double)n;
  S150:
    T3 = 15*eps;
    bgrat(&b0,&a0,&y0,&x0,w1,&T3,&ierr1);
    *w = 0.5+(0.5-*w1);
    goto End_ok;

  S160:
    n = (long)(b0);
    b0 -= (double)n;
    if(b0 == 0) {
	n -= 1;
	b0 = 1;
    }
    *w = bup(&b0,&a0,&y0,&x0,&n,&eps);
    if(x0 <= 0.7) {
	*w += bpser(&a0,&b0,&x0,&eps);
	*w1 = 0.5+(0.5-*w);
    }
    else {
	if(a0 <= 15) {
	    n = 20;
	    *w += bup(&a0,&b0,&x0,&y0,&n,&eps);
	    a0 += (double)n;
	}
	T4 = 15*eps;
	bgrat(&a0,&b0,&x0,&y0,w,&T4,&ierr1);
	*w1 = 0.5+(0.5-*w);
    }
    goto End_ok;

Use_basym:
    T5 = 100*eps;
    *w = basym (&a0,&b0,&lambda,&T5);
    *w1 = 0.5+(0.5-*w);

/*
	       TERMINATION OF THE PROCEDURE
*/
End_ok:
    if(ind != 0) { /*  swap back */
	t = *w; *w = *w1; *w1 = t;
    }
    /*
    DBGprt3(" normal --> (%g,%g)\n", *w, *w1);
    */
    return;

x_is_0:
    if(*a == 0) { *ierr = 6; return;}
b_is_0:
    *w = 0; *w1 = 1; return;

y_is_0:
    if(*b == 0) { *ierr = 7; return;}
a_is_0:
    *w = 1; *w1 = 0; return;


L_ab_small:
/*
	   PROCEDURE FOR A AND B < 1.E-3*EPS
*/
    *w = *b/(*a+*b);
    *w1 = *a/(*a+*b);
    /*
    DBGprt3(" SMALL --> (%g,%g)\n", *w, *w1);
    */
    return;
} /* bratio */

double brcmp1(int *mu,double *a,double *b,double *x,double *y)
/*
-----------------------------------------------------------------------
	  EVALUATION OF	 EXP(MU) * (X^A*Y^B / BETA(A,B))
-----------------------------------------------------------------------
*/
{
    double brcmp1,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
    int i,n;
    double T1,T2,T3,T4;

    a0 = fmin2(*a,*b);

    if(a0 < 8) {
	if(*x <= 0.375) {
	    lnx = log(*x);
	    T1 = -*x;
	    lny = alnrel(&T1);
	}
	else {
	    if(*y <= 0.375) {
		T2 = -*y;
		lnx = alnrel(&T2);
		lny = log(*y);
	    }
	    else {
		lnx = log(*x);
		lny = log(*y);
	    }
	}

	z = *a*lnx+*b*lny;
	if(a0 >= 1) {
	    z -= betaln(a,b);
	    return esum(mu,&z);
	}
	/*
	  --------------------------------------------------------------
	  PROCEDURE FOR A < 1 OR B < 1
	  --------------------------------------------------------------
	*/
	b0 = fmax2(*a,*b);

	if(b0 < 8) {

	    if(b0 <= 1) {/*		algorithm for  b0 <= 1 */

		brcmp1 = esum(mu,&z);
		if(brcmp1 == 0) return brcmp1;
		apb = *a+*b;
		if(apb <= 1) {
		    z = 1+gam1(&apb);
		}
		else {
		    u = *a+*b-1.e0;
		    z = (1+gam1(&u))/apb;
		}
		c = (1+gam1(a))*(1+gam1(b))/z;
		return brcmp1*(a0*c)/(1+a0/b0);
	    }
	    else {/*		algorithm for 1 < B0 < 8 */

		u = gamln1(&a0);
		n = (long)(b0 - 1);
		if(n >= 1) {
		    c = 1;
		    for(i=1; i<=n; i++) {
			b0 -= 1;
			c *= (b0/(a0+b0));
		    }
		    u = log(c)+u;
		}
		z -= u;
		b0 -= 1;
		apb = a0+b0;
		if(apb <= 1) {
		    t = 1+gam1(&apb);
		}
		else {
		    u = a0+b0-1.e0;
		    t = (1+gam1(&u))/apb;
		}
		return a0*esum(mu,&z)*(1+gam1(&b0))/t;
	    }
	}
	else {/*		algorithm for B0 >= 8 */

	    u = gamln1(&a0) + algdiv(&a0,&b0);
	    T3 = z-u;
	    return a0*esum(mu,&T3);
	}
    }
    else {
/*
  -----------------------------------------------------------------------
  PROCEDURE FOR A >= 8 AND B >= 8
  -----------------------------------------------------------------------
*/
	if(*a <= *b) {
	    h = *a/ *b;
	    x0 = h/(1+h);
	    y0 = 1/(1+h);
	    lambda = *a-(*a+*b)**x;
	}
	else {
	    h = *b/ *a;
	    x0 = 1/(1+h);
	    y0 = h/(1+h);
	    lambda = (*a+*b)**y-*b;
	}
	e = -(lambda/ *a);
	if(fabs(e) <= 0.6)
	    u = rlog1(&e);
	else
	    u = e-log(*x/x0);

	e = lambda/ *b;
	if(fabs(e) <= 0.6)
	    v = rlog1(&e);
	else
	    v = e-log(*y/y0);

	T4 = -(*a*u+*b*v);
	z = esum(mu,&T4);

	return M_1_SQRT_2PI*sqrt(*b*x0)*z*exp(-bcorr(a,b));
    }
}

double brcomp(double *a,double *b,double *x,double *y)
/*
-----------------------------------------------------------------------
	       EVALUATION OF  X^A * Y^B / BETA(A,B)
-----------------------------------------------------------------------
*/
{
    double brcomp,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
    int i,n;
    double T1,T2;

    if(*x == 0 || *y == 0) return 0.;

    a0 = fmin2(*a,*b);

    if(a0 < 8) {

	if(*x > 0.375) { /* goto S10 */

	    if(*y > 0.375) { /* goto S20 */
		lnx = log(*x);
		lny = log(*y);
	    }
	    else {
		T2 = -*y;
		lnx = alnrel(&T2);
		lny = log(*y);
	    }
	}
	else { /* x <= .375 */
	    lnx = log(*x);
	    T1 = -*x;
	    lny = alnrel(&T1);
	}

	z = *a*lnx+*b*lny;
	if(a0 > 1) {
	    z -= betaln(a,b);
	    brcomp = exp(z);
	    return brcomp;
	}

	/* --------------------------------------------------------------
	   a < 1 or b < 1
	   --------------------------------------------------------------
	*/
	b0 = fmax2(*a,*b);

	if(b0 < 8) {

	    if(b0 <= 1) {/*		b0 <= 1 */
		brcomp = exp(z);
		if(brcomp == 0) return brcomp;
		apb = *a+*b;
		if(apb > 1) { /* goto S50 */
		    u = *a+*b-1.e0;
		    z = (1+gam1(&u))/apb;
		}
		else {
		    z = 1+gam1(&apb);
		}

		c = (1+gam1(a))*(1+gam1(b))/z;
		return brcomp*(a0*c)/(1+a0/b0);
	    }
	    else {/*			1 < B0 < 8 */
		u = gamln1(&a0);
		n = (long)(b0 - 1);
		if(n >= 1) {
		    c = 1;
		    for(i=0; i<n; i++) {
			b0 -= 1;
			c *= (b0/(a0+b0));
		    }
		    u = log(c)+u;
		}
		z -= u;
		b0 -= 1;
		apb = a0+b0;
		if(apb <= 1) {
		    t = 1+gam1(&apb);
		}
		else {
		    u = a0+b0-1.e0;
		    t = (1+gam1(&u))/apb;
		}
		return a0*exp(z)*(1+gam1(&b0))/t;
	    }

	} /* if	 b0 = max(a,b) < 8 */


	/* ELSE			b0 = max(a,b) >= 8 */
	u = gamln1(&a0) + algdiv(&a0,&b0);
	return a0*exp(z-u);

    } /* end if	 min(a,b) < 8 */


/*
  -----------------------------------------------------------------------
  PROCEDURE FOR A >= 8 AND B >= 8
  -----------------------------------------------------------------------
*/
    if(*a <= *b) {
	h = *a/ *b;
	x0 = h/(1+h);
	y0 = 1/(1+h);
	lambda = *a-(*a+*b)**x;
    }
    else { /* a > b */
	h = *b/ *a;
	x0 = 1/(1+h);
	y0 = h/(1+h);
	lambda = (*a+*b)**y-*b;
    }

    e = -(lambda/ *a);
    if(fabs(e) <= 0.6)
	u = rlog1(&e);
    else
	u = e-log(*x/x0);

    e = lambda/ *b;
    if(fabs(e) <= 0.6)
	v = rlog1(&e);
    else
	v = e-log(*y/y0);

    z = exp(-(*a*u+*b*v));
    return M_1_SQRT_2PI*sqrt(*b*x0)*z*exp(-bcorr(a,b));
}

double bup(double *a,double *b,double *x,double *y,int *n,double *eps)
/*
-----------------------------------------------------------------------
     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
     EPS IS THE TOLERANCE USED.
-----------------------------------------------------------------------
*/
{
    double bup,ap1,apb,d,l,r,t,w;
    int k,mu,nm1;
    /*
      Obtain the scaling factor exp(-mu) and exp(mu)*(x^a*y^b/beta(a,b))/a
    */
    apb = *a+*b;
    ap1 = *a+1;
    if(*n == 1 || *a < 1 || apb < 1.1*ap1) { /* goto S10: */
	mu = 0;
	d = 1;
    }
    else {
	mu = (long)(fabs(exparg(1)));
	k  = (long)(exparg(0));
	if(k < mu) mu = k;
	t = mu;
	d = exp(-t);
    }

    bup = brcmp1(&mu,a,b,x,y)/ *a;
    if(*n == 1 || bup == 0) return bup;
    nm1 = *n-1;
    w = d;

    /*		Let k  be the index of the maximum term */
    k = 0;
    if(*b > 1) {
	if(*y <= 1.e-4) {
	    k = nm1;
	}
	else {
	    r = (*b-1)* *x/ *y - *a;
	    if(r < 1) goto S50;
	    t = nm1;
	    if(r < t) k = (long)(r); else k = (long)(t);
	}

	/*	Add the increasing terms of the series */

	for(l=0; l<k; l++) {
	    d = (apb+l)/(ap1+l)**x*d;
	    w += d;
	}
	if(k == nm1) goto L_end;
    }

 S50: /*	Add the remaining terms of the series */

    for(l=k; l<nm1; l++) {
	d = (apb+l)/(ap1+l)**x*d;
	w += d;
	if(d <= *eps*w) break;
    }

 L_end: /*	Terminate the procedure */
    bup *= w;
    return bup;
} /* bup() */

/*===============  cdf***() and cum***() etc  ===============*/

/*-- xx_check() Macros for argument checking
 *   ========== appearing  MANY times in the code : */
#define which_check_1(KK) \
    if(*which < 1 || *which > KK) {\
	if(*which < 1)\
	    *bound = 1;\
	else\
	    *bound = KK;\
	*status = -1;	return;\
    }

#define I01_check(WH,PP,STAT) \
    if(*which != WH && (PP < 0 || PP > 1)) { /* P */ \
	if(PP < 0)  *bound = 0;\
	else	    *bound = 1;\
	*status = STAT; return;\
    }

#define I001_check(WH,PP,STAT) \
    if(*which != WH && (PP <= 0 || PP > 1)) { /* P */ \
	if(PP <= 0) *bound = 0;\
	else	    *bound = 1;\
	*status = STAT; return;\
    }

#define PpQ_check \
    if(*which != 1) { /* P + Q */ \
	pq = *p+*q;\
	if(fabs(pq-0.5-0.5) > 3*spmpar(1)) {\
	    if(pq < 0)	*bound = 0;\
	    else	*bound = 1;\
	    *status = 3; return;\
	}\
    }

#define qsmall(x,S, eABS, eREL) ((S) < eABS || (x) < eREL * (S))

/*===========================================================================*/


void cdfbet(int *which,double *p,double *q,double *x,double *y,
	    double *a,double *b,int *status,double *bound)
{
/**********************************************************************

      void cdfbet(int *which,double *p,double *q,double *x,double *y,
	    double *a,double *b,int *status,double *bound)

	       Cumulative Distribution Function
			 BETa Distribution


			      Function


     Calculates any one parameter of the beta distribution given
     values for the others.


			      Arguments


     WHICH --> Integer indicating which of the next four argument
	       values is to be calculated from the others.
	       Legal range: 1..4
	       iwhich = 1 : Calculate P and Q from X,Y,A and B
	       iwhich = 2 : Calculate X and Y from P,Q,A and B
	       iwhich = 3 : Calculate A from P,Q,X,Y and B
	       iwhich = 4 : Calculate B from P,Q,X,Y and A

     P <--> The integral from 0 to X of the chi-square
	    distribution.
	    Input range: [0, 1].

     Q <--> 1-P.
	    Input range: [0, 1].
	    P + Q = 1.

     X <--> Upper limit of integration of beta density.
	    Input range: [0,1].
	    Search range: [0,1]

     Y <--> 1-X.
	    Input range: [0,1].
	    Search range: [0,1]
	    X + Y = 1.

     A <--> The first parameter of the beta density.
	    Input range: (0, +infinity).
	    Search range: [1D-100,1D100]

     B <--> The second parameter of the beta density.
	    Input range: (0, +infinity).
	    Search range: [1D-100,1D100]

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1
		4 if X + Y != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Cumulative distribution function  (P)  is calculated directly by
     code associated with the following reference.

     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
     Digit Computation of the Incomplete  Beta	Function Ratios.  ACM
     Trans. Math.  Softw. 18 (1993), 360-373.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.


			      Note


     The beta density is proportional to
	       t^(A-1) * (1-t)^(B-1)

**********************************************************************/
#define tol 1e-8
#define atol 1e-50
#define zero 1e-100
#define inf 1e100
#define one 1


static double K2 = 0;
static double K3 = 1;
static double K8 = 0.5;
static double K9 = 5;
static double fx,xhi,xlo,cum,ccum,xy,pq;
static logical qhi,qleft,qporq;
static double T4,T5,T6,T7,T10,T11,T12,T13,T14,T15;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check arguments
*/
    which_check_1(4);
    I01_check(1,*p,-2);
    I01_check(1,*q,-3);

    if(*which == 2) goto S150;
/*
     X
*/
    if(!(*x < 0 || *x > 1)) goto S140;
    if(!(*x < 0)) goto S120;
    *bound = 0;
    goto S130;
S120:
    *bound = 1;
S130:
    *status = -4;    return;

S150:
S140:
    if(*which == 2) goto S190;
/*
     Y
*/
    if(!(*y < 0 || *y > 1)) goto S180;
    if(!(*y < 0)) *bound = 1;
    *bound = 0;
    goto S170;


S170: *status = -5;    return;
S190:
S180:
    if(*which == 3) goto S210;
/*
     A
*/
    if(!(*a <= 0)) goto S200;
    *bound = 0;
    *status = -6;
    return;
S210:
S200:
    if(*which == 4) goto S230;
/*
     B
*/
    if(!(*b <= 0)) goto S220;
    *bound = 0;
    *status = -7;
    return;
S230:
S220:
    PpQ_check;

    if(*which == 2) goto S310;
/*
     X + Y
*/
    xy = *x+*y;
    if(!(fabs(xy-0.5-0.5) > 3*spmpar(1))) goto S300;
    if(!(xy < 0)) goto S280;
    *bound = 0;
    goto S290;
S280:
    *bound = 1;
S290:
    *status = 4;
    return;
S310:
S300:
    if(!(*which == 1)) qporq = *p <= *q;
/*
     Select the minimum of P or Q
     Calculate ANSWERS
*/
    if(1 == *which) {
/*
     Calculating P and Q
*/
	cumbet(x,y,a,b,p,q);
	*status = 0;
    }
    else if(2 == *which) {
/*
     Calculating X and Y
*/
	T4 = atol;
	T5 = tol;
	dstzr(&K2,&K3,&T4,&T5);
	if(!qporq) goto S340;
	*status = 0;
	dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
	*y = one-*x;
S320:
	if(!(*status == 1)) goto S330;
	cumbet(x,y,a,b,&cum,&ccum);
	fx = cum-*p;
	dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
	*y = one-*x;
	goto S320;
S330:
	goto S370;
S340:
	*status = 0;
	dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
	*x = one-*y;
S350:
	if(!(*status == 1)) goto S360;
	cumbet(x,y,a,b,&cum,&ccum);
	fx = ccum-*q;
	dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
	*x = one-*y;
	goto S350;
S370:
S360:
	if(!(*status == -1)) goto S400;
	if(!qleft) goto S380;
	*status = 1;
	*bound = 0;
	goto S390;
S380:
	*status = 2;
	*bound = 1;
S400:
S390:
	;
    }
    else if(3 == *which) {
/*
     Computing A
*/
	*a = 5;
	T6 = zero;
	T7 = inf;
	T10 = atol;
	T11 = tol;
	dstinv(&T6,&T7,&K8,&K8,&K9,&T10,&T11);
	*status = 0;
	dinvr(status,a,&fx,&qleft,&qhi);
S410:
	if(!(*status == 1)) goto S440;
	cumbet(x,y,a,b,&cum,&ccum);
	if(!qporq) goto S420;
	fx = cum-*p;
	goto S430;
S420:
	fx = ccum-*q;
S430:
	dinvr(status,a,&fx,&qleft,&qhi);
	goto S410;
S440:
	if(!(*status == -1)) goto S470;
	if(!qleft) goto S450;
	*status = 1;
	*bound = zero;
	goto S460;
S450:
	*status = 2;
	*bound = inf;
S470:
S460:
	;
    }
    else if(4 == *which) {
/*
     Computing B
*/
	*b = 5;
	T12 = zero;
	T13 = inf;
	T14 = atol;
	T15 = tol;
	dstinv(&T12,&T13,&K8,&K8,&K9,&T14,&T15);
	*status = 0;
	dinvr(status,b,&fx,&qleft,&qhi);
S480:
	if(!(*status == 1)) goto S510;
	cumbet(x,y,a,b,&cum,&ccum);
	if(!qporq) goto S490;
	fx = cum-*p;
	goto S500;
S490:
	fx = ccum-*q;
S500:
	dinvr(status,b,&fx,&qleft,&qhi);
	goto S480;
S510:
	if(!(*status == -1)) goto S540;
	if(!qleft) goto S520;
	*status = 1;
	*bound = zero;
	goto S530;
S520:
	*status = 2;
	*bound = inf;
S530:
	;
    }
S540:
    return;
#undef tol
#undef atol
#undef zero
#undef inf
#undef one
} /* cdfbet() */

void cdfbin(int *which,double *p,double *q,double *s,double *xn,
	    double *pr,double *ompr,int *status,double *bound)
{
/**********************************************************************

      void cdfbin(int *which,double *p,double *q,double *s,double *xn,
	    double *pr,double *ompr,int *status,double *bound)

	       Cumulative Distribution Function
			 BINomial distribution


			      Function


     Calculates any one parameter of the binomial
     distribution given values for the others.


			      Arguments


     WHICH --> Integer indicating which of the next four argument
	       values is to be calculated from the others.
	       Legal range: 1..4
	       iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
	       iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
	       iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
	       iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN

     P <--> The cumulation from 0 to S of the binomial distribution.
	    (Probablility of S or fewer successes in XN trials each
	    with probability of success PR.)
	    Input range: [0,1].

     Q <--> 1-P.
	    Input range: [0, 1].
	    P + Q = 1.

     S <--> The number of successes observed.
	    Input range: [0, XN]
	    Search range: [0, XN]

     XN	 <--> The number of binomial trials.
	      Input range: (0, +infinity).
	      Search range: [1E-100, 1E100]

     PR	 <--> The probability of success in each binomial trial.
	      Input range: [0,1].
	      Search range: [0,1]

     OMPR  <--> 1-PR
	      Input range: [0,1].
	      Search range: [0,1]
	      PR + OMPR = 1

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1
		4 if PR + OMPR != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula  26.5.24	 of   Abramowitz  and	 Stegun,  Handbook   of
     Mathematical   Functions (1966) is	  used	to reduce the  binomial
     distribution  to  the  cumulative incomplete    beta distribution.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.


**********************************************************************/
#define atol 1e-50
#define tol 1e-8
#define zero 1e-100
#define inf 1e100
#define one 1
static double K2 = 0;
static double K3 = 0.5;
static double K4 = 5;
static double K11 = 1;
static double fx,xhi,xlo,cum,ccum,pq,prompr;
static logical qhi,qleft,qporq;
static double T5,T6,T7,T8,T9,T10,T12,T13;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check arguments
*/
    which_check_1(4);
    I01_check(1,*p,-2);
    I01_check(1,*q,-3);


    if(*which != 3 && *xn <= 0) {/* XN */
	*bound = 0; *status = -5; return;
    }

    if((*which != 2 && *s <  0 ) ||
       (*which != 3 && *s > *xn)) {
	if(*s < 0)
	    *bound = 0;
	else
	    *bound = *xn;
	*status = -4;	 return;
    }

    I01_check(4,*pr,-6);
    I01_check(4,*ompr,-7);

    PpQ_check;

    if(*which != 4) {
	/*
	  PR + OMPR
	*/
	prompr = *pr+*ompr;
	if(fabs(prompr-0.5-0.5) > 3*spmpar(1)) {
	    if(prompr < 0)
		*bound = 0;
	    else
		*bound = 1;
	    *status = 4; return;
	}
    }

    if(*which != 1) qporq = *p <= *q;/* Select the minimum of P or Q */

/*
    Calculate ANSWERS
*/
    if(1 == *which) {
/*
     Calculating P
*/
	cumbin(s,xn,pr,ompr,p,q);
	*status = 0;
    }
    else if(2 == *which) {
/*
     Calculating S
*/
	*s = 5;
	T5 = atol;
	T6 = tol;
	dstinv(&K2,xn,&K3,&K3,&K4, &T5,&T6);
	*status = 0;
	dinvr(status,s,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumbin(s,xn,pr,ompr,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,s,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    } else {
		*status = 2; *bound = *xn;
	    }
	}
    }
    else if(3 == *which) {
/*
     Calculating XN
*/
	*xn = 5;
	T7 = zero;
	T8 = inf;
	T9 = atol;
	T10 = tol;
	dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
	*status = 0;
	dinvr(status,xn,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumbin(s,xn,pr,ompr,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,xn,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = zero;
	    } else {
		*status = 2; *bound = inf;
	    }
	}
    }
    else if(4 == *which) {
/*
     Calculating PR and OMPR
*/
	T12 = atol;
	T13 = tol;
	dstzr(&K2,&K11,&T12,&T13);
	if(qporq) {
	    *status = 0;
	    dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
	    *ompr = one-*pr;

	    while(*status == 1) {
		cumbin(s,xn,pr,ompr,&cum,&ccum);
		fx = cum-*p;
		dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
		*ompr = one-*pr;
	    }
	}
	else {

	    *status = 0;
	    dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
	    *pr = one-*ompr;

	    while(*status == 1) {
		cumbin(s,xn,pr,ompr,&cum,&ccum);
		fx = ccum-*q;
		dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
		*pr = one-*ompr;
	    }
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    } else {
		*status = 2; *bound = 1;
	    }
	}
    }
    return;
#undef atol
#undef tol
#undef zero
#undef inf
#undef one
} /* cdfbin() */

void cdfchi(int *which,double *p,double *q,double *x,double *df,
	    int *status,double *bound)
{
/**********************************************************************

      void cdfchi(int *which,double *p,double *q,double *x,double *df,
	    int *status,double *bound)

	       Cumulative Distribution Function
	       CHI-Square distribution


			      Function


     Calculates any one parameter of the chi-square
     distribution given values for the others.


			      Arguments


     WHICH --> Integer indicating which of the next three argument
	       values is to be calculated from the others.
	       Legal range: 1..3
	       iwhich = 1 : Calculate P and Q from X and DF
	       iwhich = 2 : Calculate X from P,Q and DF
	       iwhich = 3 : Calculate DF from P,Q and X

     P <--> The integral from 0 to X of the chi-square
	    distribution.
	    Input range: [0, 1].

     Q <--> 1-P.
	    Input range: (0, 1].
	    P + Q = 1.

     X <--> Upper limit of integration of the non-central
	    chi-square distribution.
	    Input range: [0, +infinity).
	    Search range: [0,1E100]

     DF <--> Degrees of freedom of the
	     chi-square distribution.
	     Input range: (0, +infinity).
	     Search range: [ 1E-100, 1E100]

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1
	       10 indicates error returned from cumgam.	 See
		  references in cdfgam

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula	26.4.19	  of Abramowitz	 and	 Stegun, Handbook  of
     Mathematical Functions   (1966) is used   to reduce the chisqure
     distribution to the incomplete distribution.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.

**********************************************************************/
#define tol 1e-8
#define atol 1e-50
#define zero 1e-100
#define inf 1e100
static double K2 = 0;
static double K4 = 0.5;
static double K5 = 5;
static double fx,cum,ccum,pq,porq;
static logical qhi,qleft,qporq;
static double T3,T6,T7,T8,T9,T10,T11;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check arguments
*/
    which_check_1(3);
    I01_check(1,*p,-2);
    I01_check(1,*q,-3);

    if(*which != 2 && *x < 0) { /* X */
	*bound = 0;
	*status = -4; return;
    }
    if(*which == 3 && *df <= 0) { /* DF */
	*bound = 0;
	*status = -5; return;
    }

    PpQ_check;

    if(*which != 1) {
	/* Select the minimum of P or Q */
	qporq = *p <= *q;
	if(qporq)
	    porq = *p;
	else
	    porq = *q;
    }

/* ------------------------------------------------
     Calculate ANSWERS
*/
    switch(*which) {
    case 1:
	/*
	  Calculating P and Q
	*/
	*status = 0;
	cumchi(x,df,p,q);
	if(porq > 1.5) {
	    *status = 10;
	    return;
	}
	break;

    case 2:
	/*
	  Calculating X
	*/
	*x = 5;
	T3 = inf;
	T6 = atol;
	T7 = tol;
	dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
	*status = 0;
	dinvr(status,x,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumchi(x,df,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    if(fx+porq > 1.5) {
		*status = 10; return;
	    }
	    dinvr(status,x,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    } else {
		*status = 2; *bound = inf;
	    }
	}

	break;

    case 3:
	/*
	  Calculating DF
	*/
	*df = 5;
	T8 = zero;
	T9 = inf;
	T10 = atol;
	T11 = tol;
	dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
	*status = 0;
	dinvr(status,df,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumchi(x,df,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;

	    if(fx+porq > 1.5) {
		*status = 10; return;
	    }
	    dinvr(status,df,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = zero;
	    } else {
		*status = 2; *bound = inf;
	    }
	}
    }
    return;
#undef tol
#undef atol
#undef zero
#undef inf
} /* cdfchi() */

void cdfchn(int *which,double *p,double *q,double *x,double *df,
	    double *pnonc,int *status,double *bound)
{
/**********************************************************************

      void cdfchn(int *which,double *p,double *q,double *x,double *df,
	    double *pnonc,int *status,double *bound)

	       Cumulative Distribution Function
	       Non-central Chi-Square


			      Function


     Calculates any one parameter of the non-central chi-square
     distribution given values for the others.


			      Arguments


     WHICH --> Integer indicating which of the next three argument
	       values is to be calculated from the others.
	       Input range: 1..4
	       iwhich = 1 : Calculate P and Q from X and DF
	       iwhich = 2 : Calculate X from P,DF and PNONC
	       iwhich = 3 : Calculate DF from P,X and PNONC
	       iwhich = 3 : Calculate PNONC from P,X and DF

     P <--> The integral from 0 to X of the non-central chi-square
	    distribution.
	    Input range: [0, 1-1E-16).

     Q <--> 1-P.
	    Q is not used by this subroutine and is only included
	    for similarity with other cdf* routines.

     X <--> Upper limit of integration of the non-central
	    chi-square distribution.
	    Input range: [0, +infinity).
	    Search range: [0,1E100]

     DF <--> Degrees of freedom of the non-central
	     chi-square distribution.
	     Input range: (0, +infinity).
	     Search range: [ 1E-100, 1E100]

     PNONC <--> Non-centrality parameter of the non-central
		chi-square distribution.
		Input range: [0, +infinity).
		Search range: [0,1E4]

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula  26.4.25	of   Abramowitz	  and	Stegun,	 Handbook  of
     Mathematical  Functions (1966) is used to compute the cumulative
     distribution function.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.


			    WARNING

     The computation time  required for this  routine is proportional
     to the noncentrality  parameter  (PNONC).	Very large  values of
     this parameter can consume immense	 computer resources.  This is
     why the search range is bounded by 10,000.

**********************************************************************/
#define tent4 1e4
#define tol 1e-8
#define atol 1e-50
#define zero 1e-100
#define one ( 1 - 1e-16 )
#define inf 1e100
static double K1 = 0;
static double K3 = 0.5;
static double K4 = 5;
static double fx,cum,ccum;
static logical qhi,qleft;
static double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check arguments
*/

    which_check_1(4);
    I01_check(1,*p,-2);
    I01_check(1,*q,-3);

    if(*which != 2 && *x < 0) {
	*bound = 0; *status = -4; return;
    }
    if(*which != 3 && *df <= 0) {
	*bound = 0; *status = -5; return;
    }
    if(*which != 4 && *pnonc < 0) {
	*bound = 0; *status = -6; return;
    }
/*
     Calculate ANSWERS
*/
    if(1 == *which) {
/*
     Calculating P and Q
*/
	cumchn(x,df,pnonc,p,q);
	*status = 0;
    }
    else if(2 == *which) {
/*
     Calculating X
*/
	*x = 5;
	T2 = inf;
	T5 = atol;
	T6 = tol;
	dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
	*status = 0;
	dinvr(status,x,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumchn(x,df,pnonc,&cum,&ccum);
	    fx = cum-*p;
	    dinvr(status,x,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) { *status = 1; *bound = 0;
	    }
	    else {	*status = 2; *bound = inf;}
	}

    }
    else if(3 == *which) {
/*
     Calculating DF
*/
	*df = 5;
	T7 = zero;
	T8 = inf;
	T9 = atol;
	T10 = tol;
	dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
	*status = 0;
	dinvr(status,df,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumchn(x,df,pnonc,&cum,&ccum);
	    fx = cum-*p;
	    dinvr(status,df,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = zero;
	    } else {
		*status = 2; *bound = inf;
	    }
	}
    }
    else if(4 == *which) {
/*
     Calculating PNONC
*/
	*pnonc = 5;
	T11 = tent4;
	T12 = atol;
	T13 = tol;
	dstinv(&K1,&T11,&K3,&K3,&K4,&T12,&T13);
	*status = 0;
	dinvr(status,pnonc,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumchn(x,df,pnonc,&cum,&ccum);
	    fx = cum-*p;
	    dinvr(status,pnonc,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = zero;
	    } else {
		*status = 2; *bound = tent4;
	    }
	}
    }
    return;
#undef tent4
#undef tol
#undef atol
#undef zero
#undef one
#undef inf
} /* cdfchn() */

void cdff(int *which,double *p,double *q,double *f,double *dfn,
	  double *dfd,int *status,double *bound)
{
/**********************************************************************

      void cdff(int *which,double *p,double *q,double *f,double *dfn,
		double *dfd,int *status,double *bound)

		   Cumulative Distribution Function
			    F distribution

			      Function


     Calculates any one parameter of the F distribution
     given values for the others.


			      Arguments


     WHICH --> Integer indicating which of the next four argument
	       values is to be calculated from the others.
	       Legal range: 1..4
	       iwhich = 1 : Calculate P and Q from F,DFN and DFD
	       iwhich = 2 : Calculate F from P,Q,DFN and DFD
	       iwhich = 3 : Calculate DFN from P,Q,F and DFD
	       iwhich = 4 : Calculate DFD from P,Q,F and DFN

       P <--> The integral from 0 to F of the f-density.
	      Input range: [0,1].

       Q <--> 1-P.
	      Input range: (0, 1].
	      P + Q = 1.

       F <--> Upper limit of integration of the f-density.
	      Input range: [0, +infinity).
	      Search range: [0,1E100]

     DFN < --> Degrees of freedom of the numerator sum of squares.
	       Input range: (0, +infinity).
	       Search range: [ 1E-100, 1E100]

     DFD < --> Degrees of freedom of the denominator sum of squares.
	       Input range: (0, +infinity).
	       Search range: [ 1E-100, 1E100]

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula   26.6.2	of   Abramowitz	  and	Stegun,	 Handbook  of
     Mathematical  Functions (1966) is used to reduce the computation
     of the  cumulative	 distribution function for the	F  variate to
     that of an incomplete beta.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.

			      WARNING

     The value of the  cumulative  F distribution is  not necessarily
     monotone in  either degrees of freedom.  There  thus may  be two
     values  that  provide a given CDF	value.	 This routine assumes
     monotonicity and will find an arbitrary one of the two values.

**********************************************************************/
#define tol 1e-8
#define atol 1e-50
#define zero 1e-100
#define inf 1e100

static double K0 = 0;
static double Half = 0.5;
static double K5 = 5;

    double pq,fx,cum,ccum;
    logical qhi,qleft,qporq=0;
    double T3,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15;

/*
     Check arguments
*/
    which_check_1(4);
    I01_check(1,*p,-2);
    I001_check(1,*q,-3);

    if(*which != 2 && *f < 0) {
	*bound = 0; *status = -4; return;
    }
    if(*which != 3 && *dfn <= 0) {
	*bound = 0; *status = -5; return;
    }
    if(*which != 4 && *dfd <= 0) {
	*bound = 0; *status = -6; return;
    }

    PpQ_check;

    if(*which != 1) qporq = *p <= *q;
/*
     Select the minimum of P or Q

     Calculate ANSWERS
*/

    switch(*which) {
    case 1:
/*
     Calculating P
*/
	cumf(f,dfn,dfd,p,q);
	*status = 0;
	return;

    case 2:
/*
     Calculating F
*/
	*f = 5;
	T3 = inf;
	T6 = atol;
	T7 = tol;
	dstinv(&K0,&T3,&Half,&Half,&K5,&T6,&T7);
	*status = 0;
	dinvr(status,f,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumf(f,dfn,dfd,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,f,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;


    case 3:
/*
     Calculating DFN
*/
	*dfn = 5;
	T8 = zero;
	T9 = inf;
	T10 = atol;
	T11 = tol;
	dstinv(&T8,&T9,&Half,&Half,&K5,&T10,&T11);
	*status = 0;
	dinvr(status,dfn,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumf(f,dfn,dfd,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,dfn,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = zero;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;


    case 4:
/*
     Calculating DFD
*/
	*dfd = 5;
	T12 = zero;
	T13 = inf;
	T14 = atol;
	T15 = tol;
	dstinv(&T12,&T13,&Half,&Half,&K5,&T14,&T15);
	*status = 0;
	dinvr(status,dfd,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumf(f,dfn,dfd,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,dfn,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = zero;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;
    }

#undef tol
#undef atol
#undef zero
#undef inf
} /* cdff() */

void cdffnc(int *which,double *p,double *q,double *f,double *dfn,
	    double *dfd,double *pnonc,int *status,double *bound)
{
/**********************************************************************

      void cdffnc(int *which,double *p,double *q,double *f,double *dfn,
	    double *dfd,double *pnonc,int *status,double *bound)

	       Cumulative Distribution Function
	       Non-central F distribution


			      Function


     Calculates any one parameter of the Non-central F
     distribution given values for the others.


			      Arguments


     WHICH --> Integer indicating which of the next five argument
	       values is to be calculated from the others.
	       Legal range: 1..5
	       iwhich = 1 : Calculate P and Q from F,DFN,DFD and PNONC
	       iwhich = 2 : Calculate F from P,Q,DFN,DFD and PNONC
	       iwhich = 3 : Calculate DFN from P,Q,F,DFD and PNONC
	       iwhich = 4 : Calculate DFD from P,Q,F,DFN and PNONC
	       iwhich = 5 : Calculate PNONC from P,Q,F,DFN and DFD

       P <--> The integral from 0 to F of the non-central f-density.
	      Input range: [0,1-1E-16).

       Q <--> 1-P.
	      Q is not used by this subroutine and is only included
	      for similarity with other cdf* routines.

       F <--> Upper limit of integration of the non-central f-density.
	      Input range: [0, +infinity).
	      Search range: [0,1E100]

     DFN < --> Degrees of freedom of the numerator sum of squares.
	       Input range: (0, +infinity).
	       Search range: [ 1E-100, 1E100]

     DFD < --> Degrees of freedom of the denominator sum of squares.
	       Must be in range: (0, +infinity).
	       Input range: (0, +infinity).
	       Search range: [ 1E-100, 1E100]

     PNONC <-> The non-centrality parameter
	       Input range: [0,infinity)
	       Search range: [0,1E4]

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula  26.6.20	of   Abramowitz	  and	Stegun,	 Handbook  of
     Mathematical  Functions (1966) is used to compute the cumulative
     distribution function.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.

			    WARNING

     The computation time  required for this  routine is proportional
     to the noncentrality  parameter  (PNONC).	Very large  values of
     this parameter can consume immense	 computer resources.  This is
     why the search range is bounded by 10,000.

			      WARNING

     The  value	 of the	 cumulative  noncentral F distribution is not
     necessarily monotone in either degrees  of freedom.  There	 thus
     may be two values that provide a given  CDF value.	 This routine
     assumes monotonicity  and will find  an arbitrary one of the two
     values.

**********************************************************************/
#define tent4 1e4
#define tol 1e-8
#define atol 1e-50
#define Azero 1e-100
#define inf 1e100
static double K0 = 0;
static double Half = 0.5;
static double K5 = 5;

    double fx,cum,ccum;
    logical qhi,qleft;
    double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17;

/*
     Check arguments
*/
    which_check_1(5);
    I01_check(1,*p,-2);
    I01_check(1,*q,-3);

    if(*which != 2 && *f < 0) {
	*bound = 0; *status = -4; return;
    }
    if(*which != 3 && *dfn <= 0) {
	*bound = 0; *status = -5; return;
    }
    if(*which != 4 && *dfd <= 0) {
	*bound = 0; *status = -6; return;
    }
    if(*which != 5 && *pnonc < 0) {
	*bound = 0; *status = -7; return;
    }

/*
     Calculate ANSWERS
*/
    switch(*which) {
    case 1:
/*
  Calculating P
*/
	cumfnc(f,dfn,dfd,pnonc,p,q);
	*status = 0;
	return;

    case 2:
/*
  Calculating F
*/
	*f = 5;
	T2 = inf;
	T5 = atol;
	T6 = tol;
	dstinv(&K0,&T2,&Half,&Half,&K5,&T5,&T6);
	*status = 0;
	dinvr(status,f,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
	    fx = cum-*p;
	    dinvr(status,f,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;

    case 3:
/*
  Calculating DFN
*/
	*dfn = 5;
	T7 = Azero;
	T8 = inf;
	T9 = atol;
	T10 = tol;
	dstinv(&T7,&T8,&Half,&Half,&K5,&T9,&T10);
	*status = 0;
	dinvr(status,dfn,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
	    fx = cum-*p;
	    dinvr(status,dfn,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = Azero;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;

    case 4:
/*
  Calculating DFD
*/
	*dfd = 5;
	T11 = Azero;
	T12 = inf;
	T13 = atol;
	T14 = tol;
	dstinv(&T11,&T12,&Half,&Half,&K5,&T13,&T14);
	*status = 0;
	dinvr(status,dfd,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
	    fx = cum-*p;
	    dinvr(status,dfd,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = Azero;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;

    case 5:
/*
  Calculating PNONC
*/
	*pnonc = 5;
	T15 = tent4;
	T16 = atol;
	T17 = tol;
	dstinv(&K0,&T15,&Half,&Half,&K5,&T16,&T17);
	*status = 0;
	dinvr(status,pnonc,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
	    fx = cum-*p;
	    dinvr(status,pnonc,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    }
	    else {
		*status = 2; *bound = tent4;
	    }
	}
	return;
    }
#undef tent4
#undef tol
#undef atol
#undef Azero
#undef inf
} /* cdffnc() */

void cdfgam(int *which, double *p, double *q,
	    double *x, double *shape, double *scale,
	    int *status, double *bound)
{
/**********************************************************************

		Cumulative Distribution Function
			GAMma Distribution


			Function


	Calculates any one parameter of the gamma
	distribution given values for the others.


			Arguments


    WHICH -->	Integer indicating which of the next four argument
		values is to be calculated from the others.
		Legal range: 1..4
		which = 1 : Calculate P and Q from X, SHAPE and SCALE
		which = 2 : Calculate X from P, Q, SHAPE and SCALE
		which = 3 : Calculate SHAPE from P,Q,X and SCALE
		which = 4 : Calculate SCALE from P,Q,X and SHAPE

     P	<-->	The integral from 0 to X of the gamma density.
		Input range: [0,1].

     Q	<-->	= 1-P.
		Input range: (0, 1].	P + Q = 1.

     X	<-->	The upper limit of integration of the gamma density.
		Input range: [0, +infinity).
		Search range: [0,1E100]

     SHAPE <--> The shape parameter of the gamma density.
		Input range: (0, +infinity).
		Search range: [1E-100,1E100]

     SCALE <--> The scale parameter of the gamma density.
		Input range: (0, +infinity).
		Search range: (1E-100,1E100]

 === CHANGED by Martin Maechler ===
 === is now the `real' scale, i.e. new scale = 1 / {old scale}


     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1
		10 if the gamma or inverse gamma routine cannot
		   compute the answer.	Usually happens only for
		   X and SHAPE very large (gt 1E10 or more)

     BOUND <--	 Undefined if STATUS is 0

		 Bound exceeded by parameter number I if STATUS
		 is negative.

		 Lower search bound if STATUS is 1.

		 Upper search bound if STATUS is 2.


			Method


	Cumulative distribution function (P) is calculated directly by
	the code associated with:

	DiDinato, A. R. and Morris, A. H. Computation of the  incomplete
	gamma function	ratios	and their	inverse.   ACM	Trans.	Math.
	Softw. 12 (1986), 377-393.

	Computation of other parameters involve a seach for a value that
	produces  the desired  value  of P.   The search relies	 on  the
	monotinicity of P with the other parameter.


			Note


	The gamma density is proportional to

	  t ^ (shape - 1) * EXP(- t / scale)

	  (meaning of `scale' changed by Martin Maechler, see above)
 **************************************************************************
*/
#define tol 1e-8
#define atol 1e-50
#define Azero 1e-100
#define inf 1e100

static double Half = 0.5;
static double K5 = 5;
static double Km1 = -1;

    double xx,fx,xscale,cum,ccum,pq,porq=0.;
    int ierr;
    logical qhi,qleft,qporq=0;
    double T3,T4,T7,T8;
/*
     Check arguments
*/
    which_check_1(4);
    I01_check (1,*p,-2);
    I001_check(1,*q,-3);

    if(*which != 2 &&  *x    <	0) { *bound = 0; *status = -4; return; }
    if(*which != 3 && *shape <= 0) { *bound = 0; *status = -5; return; }
    if(*which != 4 && *scale <= 0) { *bound = 0; *status = -6; return; }

    PpQ_check;

    if(*which != 1) { /* Select the minimum of P or Q */
	qporq = *p <= *q;
	if(qporq)
	    porq = *p;
	else
	    porq = *q;
    }

/*
     Calculate ANSWERS
*/
    switch(*which) {

    case 1:/*		Calculating P */
	*status = 0;
	xscale = *x / *scale;
	cumgam(&xscale,shape,p,q);
	/* was " if(porq > 1.5) " :
	   BUG in dcdflib-c-1.1 : porq is NOT initialized for which=1 */
	if(*p > 1.5 || *q > 1.5)
	    *status = 10;
	return;

    case 2:/*		Computing X */
	gaminv(shape,&xx,&Km1, p,q,&ierr);
	if(ierr < 0) {
	    *status = 10;
	}
	else  {
	    *x = xx * *scale;
	    *status = 0;
	}
	return;


    case 3:/*		Computing SHAPE */
	*shape = 5;
	xscale = *x / *scale;
	T3 = Azero;
	T4 = inf;
	T7 = atol;
	T8 = tol;
	dstinv(&T3,&T4,&Half,&Half,&K5,&T7,&T8);
	*status = 0;
	dinvr(status,shape,&fx,&qleft,&qhi);

	while(*status == 1) {/* inversion iterations */
	    cumgam(&xscale,shape,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;

	    if(( qporq && cum  > 1.5) ||
	       (!qporq && ccum > 1.5)) {
		*status = 10;
		return;
	    }
	    dinvr(status,shape,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) { *status = 1; *bound = Azero;
	    } else {	*status = 2; *bound = inf; }
	}
	return;


    case 4:/*		Computing SCALE */
	gaminv(shape,&xx, &Km1, p,q,&ierr);
	if(ierr < 0) {
	    *status = 10;
	}
	else  {
	    *scale = *x / xx;
	    *status = 0;
	}
	return;
    }
#undef tol
#undef atol
#undef Azero
#undef inf
} /* cdfgam() */

void cdfnbn(int *which,double *p,double *q,double *s,double *xn,
	    double *pr,double *ompr,int *status,double *bound)
{
/**********************************************************************

      void cdfnbn(int *which,double *p,double *q,double *s,double *xn,
	    double *pr,double *ompr,int *status,double *bound)

	       Cumulative Distribution Function
	       Negative BiNomial distribution


			      Function


     Calculates any one parameter of the negative binomial
     distribution given values for the others.

     The  cumulative  negative	 binomial  distribution	 returns  the
     probability that there  will be  F or fewer failures before  the
     XNth success in binomial trials each of which has probability of
     success PR.

     The individual term of the negative binomial is the probability of
     S failures before XN successes and is
	  Choose( S, XN+S-1 ) * PR^(XN) * (1-PR)^S


			      Arguments


     WHICH --> Integer indicating which of the next four argument
	       values is to be calculated from the others.
	       Legal range: 1..4
	       iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
	       iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
	       iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
	       iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN

     P <--> The cumulation from 0 to S of the  negative
	    binomial distribution.
	    Input range: [0,1].

     Q <--> 1-P.
	    Input range: (0, 1].
	    P + Q = 1.

     S <--> The upper limit of cumulation of the binomial distribution.
	    There are F or fewer failures before the XNth success.
	    Input range: [0, +infinity).
	    Search range: [0, 1E100]

     XN	 <--> The number of successes.
	      Input range: [0, +infinity).
	      Search range: [0, 1E100]

     PR	 <--> The probability of success in each binomial trial.
	      Input range: [0,1].
	      Search range: [0,1].

     OMPR  <--> 1-PR
	      Input range: [0,1].
	      Search range: [0,1]
	      PR + OMPR = 1

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1
		4 if PR + OMPR != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula   26.5.26	 of   Abramowitz  and  Stegun,	Handbook   of
     Mathematical Functions (1966) is used  to	reduce calculation of
     the cumulative distribution  function to that of  an  incomplete
     beta.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.

**********************************************************************/
#define tol 1e-8
#define atol 1e-50
#define inf 1e100
#define one 1

static double K0 = 0;
static double K1 = 1;
static double K5 = 5;
static double Half = 0.5;

    double fx,xhi,xlo,pq,prompr,cum,ccum;
    logical qhi,qleft,qporq=0;
    double T3,T6,T7,T8,T9,T10,T12,T13;
/*
     Check arguments
*/
    which_check_1(4);
    I01_check(1,*p,-2);
    I001_check(1,*q,-3);

    if(*which != 2 && *s  < 0) { *bound = 0; *status = -4; return; }
    if(*which != 3 && *xn < 0) { *bound = 0; *status = -5; return; }

    I01_check(4,*pr,-6);
    I01_check(4,*ompr,-7);

    PpQ_check;

    if(*which != 4) {
/*
  PR + OMPR
*/
	prompr = *pr+*ompr;
	if(fabs(prompr-0.5-0.5) > 3*spmpar(1)) {
	    if(prompr < 0)
		*bound = 0;
	    else
		*bound = 1;
	    *status = 4;
	    return;
	}
    }

    if(*which != 1) qporq = *p <= *q;
/*
     Select the minimum of P or Q

     Calculate ANSWERS
*/

    switch(*which) {

    case 1:
/*
  Calculating P
*/
	cumnbn(s,xn,pr,ompr,p,q);
	*status = 0;
	return;

    case 2:
/*
  Calculating S
*/
	*s = 5;
	T3 = inf;
	T6 = atol;
	T7 = tol;
	dstinv(&K0,&T3,&Half,&Half,&K5,&T6,&T7);
	*status = 0;
	dinvr(status,s,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumnbn(s,xn,pr,ompr,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,s,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;


    case 3:
/*
  Calculating XN
*/
	*xn = 5;
	T8 = inf;
	T9 = atol;
	T10 = tol;
	dstinv(&K0,&T8,&Half,&Half,&K5,&T9,&T10);
	*status = 0;
	dinvr(status,xn,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumnbn(s,xn,pr,ompr,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,s,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    }
	    else {
		*status = 2; *bound = inf;
	    }
	}
	return;


    case 4:
/*
  Calculating PR and OMPR
*/
	T12 = atol;
	T13 = tol;
	dstzr(&K0,&K1,&T12,&T13);
	*status = 0;
	if(qporq) {
	    dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
	    *ompr = one-*pr;

	    while(*status == 1) {

		cumnbn(s,xn,pr,ompr,&cum,&ccum);
		fx = cum-*p;
		dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
		*ompr = one-*pr;
	    }
	}
	else {
	    dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
	    *pr = one-*ompr;

	    while(*status == 1) {
		cumnbn(s,xn,pr,ompr,&cum,&ccum);
		fx = ccum-*q;
		dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
		*pr = one-*ompr;
	    }
	}

	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    }
	    else {
		*status = 2; *bound = 1;
	    }
	}
	return;
    }
#undef tol
#undef atol
#undef inf
#undef one
} /* cdfnbn() */

void cdfnor(int *which,double *p,double *q,double *x,double *mean,
	    double *sd,int *status,double *bound)
{
/**********************************************************************

      void cdfnor(int *which,double *p,double *q,double *x,double *mean,
	    double *sd,int *status,double *bound)

	       Cumulative Distribution Function
	       NORmal distribution


			      Function


     Calculates any one parameter of the normal
     distribution given values for the others.


			      Arguments


     WHICH  --> Integer indicating  which of the  next	parameter
     values is to be calculated using values  of the others.
     Legal range: 1..4
	       iwhich = 1 : Calculate P and Q from X,MEAN and SD
	       iwhich = 2 : Calculate X from P,Q,MEAN and SD
	       iwhich = 3 : Calculate MEAN from P,Q,X and SD
	       iwhich = 4 : Calculate SD from P,Q,X and MEAN

     P <--> The integral from -infinity to X of the normal density.
	    Input range: (0,1].

     Q <--> 1-P.
	    Input range: (0, 1].
	    P + Q = 1.

     X < --> Upper limit of integration of the normal-density.
	     Input range: ( -infinity, +infinity)

     MEAN <--> The mean of the normal density.
	       Input range: (-infinity, +infinity)

     SD <--> Standard Deviation of the normal density.
	     Input range: (0, +infinity).

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method




     A slightly modified version of ANORM from

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     is used to calulate the  cumulative standard normal distribution.

     The rational functions from pages	90-95  of Kennedy and Gentle,
     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
     starting values to Newton's Iterations which compute the inverse
     standard normal.  Therefore no  searches  are necessary for  any
     parameter.

     For X < -15, the asymptotic expansion for the normal is used  as
     the starting value in finding the inverse standard normal.
     This is formula 26.2.12 of Abramowitz and Stegun.


			      Note


      The normal density is proportional to
      exp( - 0.5 * (( X - MEAN)/SD)**2)

**********************************************************************/
    double z,pq;
/*
     Check arguments
*/
    *status = 0;

    which_check_1(4);
    I001_check(1,*p,-2);
    I001_check(1,*q,-3);

    if(*which != 4 && *sd <= 0) { *bound = 0; *status = -6; return; }

    PpQ_check ;

    switch(*which) {

    case 1: /*	   Computing P */
	z = (*x-*mean)/ *sd;
	cumnor(&z,p,q);
	return;

    case 2: /*	   Computing X */
	z = dinvnr(p,q);
	*x = *sd*z+*mean;
	return;

    case 3: /*	   Computing the MEAN */
	z = dinvnr(p,q);
	*mean = *x-*sd*z;
	return;

    case 4: /*	   Computing SD */
	z = dinvnr(p,q);
	*sd = (*x-*mean)/z;
	return;
    }
} /* cdfnor() */

void cdfpoi(int *which,double *p,double *q,double *s,double *xlam,
	    int *status,double *bound)
{
/**********************************************************************

      void cdfpoi(int *which,double *p,double *q,double *s,double *xlam,
	    int *status,double *bound)

	       Cumulative Distribution Function
	       POIsson distribution


			      Function


     Calculates any one parameter of the Poisson
     distribution given values for the others.


			      Arguments


     WHICH --> Integer indicating which	 argument
	       value is to be calculated from the others.
	       Legal range: 1..3
	       iwhich = 1 : Calculate P and Q from S and XLAM
	       iwhich = 2 : Calculate A from P,Q and XLAM
	       iwhich = 3 : Calculate XLAM from P,Q and S

	P <--> The cumulation from 0 to S of the poisson density.
	       Input range: [0,1].

	Q <--> 1-P.
	       Input range: (0, 1].
	       P + Q = 1.

	S <--> Upper limit of cumulation of the Poisson.
	       Input range: [0, +infinity).
	       Search range: [0,1E100]

     XLAM <--> Mean of the Poisson distribution.
	       Input range: [0, +infinity).
	       Search range: [0,1E100]

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula   26.4.21	of   Abramowitz	 and   Stegun,	 Handbook  of
     Mathematical Functions (1966) is used  to reduce the computation
     of	 the cumulative distribution function to that  of computing a
     chi-square, hence an incomplete gamma function.

     Cumulative	 distribution function	(P) is	calculated  directly.
     Computation of other parameters involve a seach for a value that
     produces  the desired value of  P.	  The  search relies  on  the
     monotinicity of P with the other parameter.

**********************************************************************/
#define tol 1e-8
#define atol 1e-50
#define inf 1e100

static double K0 = 0;
static double Half = 0.5;
static double K5 = 5;

    double fx,cum,ccum,pq;
    logical qhi,qleft,qporq=0;
    double T3,T6,T7,T8,T9,T10;
/*
     Check arguments
*/
    which_check_1(3);
    I01_check(1,*p,-2);
    I001_check(1,*q,-3);

    if(*which != 2 &&  *s    <	0) { *bound = 0; *status = -4; return; }
    if(*which != 3 && *xlam  <	0) { *bound = 0; *status = -5; return; }

    PpQ_check;

    if(!(*which == 1)) qporq = *p <= *q;
/*
     Select the minimum of P or Q

     Calculate ANSWERS
*/
    switch(*which) {

    case 1: /*		Calculating P */
	cumpoi(s,xlam,p,q);
	*status = 0; return;

    case 2: /*		Calculating S */
	*s = 5;
	T3 = inf;
	T6 = atol;
	T7 = tol;
	dstinv(&K0,&T3,&Half,&Half,&K5,&T6,&T7);
	*status = 0;
	dinvr(status,s,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumpoi(s,xlam,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,s,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {	*status = 1; *bound = 0; }
	    else {	*status = 2; *bound = inf; }
	}
	return;

    case 3: /*		Calculating XLAM */
	*xlam = 5;
	T8 = inf;
	T9 = atol;
	T10 = tol;
	dstinv(&K0,&T8,&Half,&Half,&K5,&T9,&T10);
	*status = 0;
	dinvr(status,xlam,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumpoi(s,xlam,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,xlam,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {	*status = 1; *bound = 0; }
	    else {	*status = 2; *bound = inf; }
	}
	return;
    }
#undef tol
#undef atol
#undef inf
} /* cdfpoi() */

void cdft(int *which,double *p,double *q,double *t,double *df,
	  int *status,double *bound)
{
/**********************************************************************

      void cdft(int *which,double *p,double *q,double *t,double *df,
	  int *status,double *bound)

	       Cumulative Distribution Function
			 T distribution


			      Function


     Calculates any one parameter of the t distribution given
     values for the others.


			      Arguments


     WHICH --> Integer indicating which	 argument
	       values is to be calculated from the others.
	       Legal range: 1..3
	       iwhich = 1 : Calculate P and Q from T and DF
	       iwhich = 2 : Calculate T from P,Q and DF
	       iwhich = 3 : Calculate DF from P,Q and T

	P <--> The integral from -infinity to t of the t-density.
	       Input range: (0,1].

	Q <--> 1-P.
	       Input range: (0, 1].
	       P + Q = 1.

	T <--> Upper limit of integration of the t-density.
	       Input range: ( -infinity, +infinity).
	       Search range: [ -1E100, 1E100 ]

	DF <--> Degrees of freedom of the t-distribution.
		Input range: (0 , +infinity).
		Search range: [1e-100, 1E10]

     STATUS <-- 0 if calculation completed correctly
	       -I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		  search bound
		2 if answer appears to be higher than greatest
		  search bound
		3 if P + Q != 1

     BOUND <-- Undefined if STATUS is 0

	       Bound exceeded by parameter number I if STATUS
	       is negative.

	       Lower search bound if STATUS is 1.

	       Upper search bound if STATUS is 2.


			      Method


     Formula  26.5.27  of   Abramowitz	 and  Stegun,	Handbook   of
     Mathematical Functions  (1966) is used to reduce the computation
     of the cumulative distribution function to that of an incomplete
     beta.

     Computation of other parameters involve a seach for a value that
     produces  the desired  value  of P.   The search relies  on  the
     monotinicity of P with the other parameter.

**********************************************************************/
#define tol 1e-8
#define atol 1e-50
#define Azero 1e-100
#define inf 1e100
#define rtinf 1e100
#define maxdf 1e10

static double Half = 0.5;
static double K5 = 5;

    double fx,cum,ccum,pq;
    logical qhi,qleft,qporq=0;
    double T2,T3,T6,T7,T8,T9,T10,T11;
/*
     Check arguments
*/
    which_check_1(3);
    I001_check(1,*p,-2);
    I001_check(1,*q,-3);

    if(*which != 3 &&  *df  <=	0) { *bound = 0; *status = -5; return; }

    PpQ_check;

    if(!(*which == 1)) qporq = *p <= *q;
/*
     Select the minimum of P or Q

     Calculate ANSWERS
*/
    switch(*which) {

    case 1: /*		Computing P and Q */
	cumt(t,df,p,q);
	*status = 0; return;


    case 2: /*		Computing T */

	/* .. Get initial approximation for T */
	*t = dt1(p,q,df);
	T2 = -rtinf;
	T3 = rtinf;
	T6 = atol;
	T7 = tol;
	dstinv(&T2,&T3,&Half,&Half,&K5,&T6,&T7);
	*status = 0;
	dinvr(status,t,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumt(t,df,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,t,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {	*status = 1; *bound = -rtinf; }
	    else {	*status = 2; *bound =  rtinf; }
	}
	return;


    case 3:/*		Computing DF */
	*df = 5;
	T8 = Azero;
	T9 = maxdf;
	T10 = atol;
	T11 = tol;
	dstinv(&T8,&T9,&Half,&Half,&K5,&T10,&T11);
	*status = 0;
	dinvr(status,df,&fx,&qleft,&qhi);
	while(*status == 1) {
	    cumt(t,df,&cum,&ccum);
	    if(qporq)
		fx = cum-*p;
	    else
		fx = ccum-*q;
	    dinvr(status,df,&fx,&qleft,&qhi);
	}

	if(*status == -1) {
	    if(qleft) {	*status = 1; *bound = Azero; }
	    else {	*status = 2; *bound = maxdf; }
	}
	return;
    }
#undef tol
#undef atol
#undef Azero
#undef inf
#undef rtinf
#undef maxdf
} /* cdft() */

void cdftnc(int *which,double *p,double *q,double *t,double *df,
	    double *pnonc,int *status,double *bound)
{
/**********************************************************************

   void cdftnc(int *which,double *p,double *q,double *t,double *df,
	       double *pnonc,int *status,double *bound)

		Cumulative Distribution Function
		   Non-Central T distribution

				Function

      Calculates any one parameter of the noncentral t distribution give
      values for the others.

				Arguments

      WHICH --> Integer indicating which  argument
		values is to be calculated from the others.
		Legal range: 1..4
		iwhich = 1 : Calculate P and Q from T,DF,PNONC
		iwhich = 2 : Calculate T from P,Q,DF,PNONC
		iwhich = 3 : Calculate DF from P,Q,T
		iwhich = 4 : Calculate PNONC from P,Q,DF,T

	 P <--> The integral from -infinity to t of the noncentral t-den
	       Input range: (0,1].

	 Q <--> 1-P.
	       Input range: (0, 1].
		P + Q = 1.

	 T <--> Upper limit of integration of the noncentral t-density.
		Input range: ( -infinity, +infinity).
		Search range: [ -1E100, 1E100 ]

	 DF <--> Degrees of freedom of the noncentral t-distribution.
		 Input range: (0 , +infinity).
		 Search range: [1e-100, 1E10]

      PNONC <--> Noncentrality parameter of the noncentral t-distributio
		 Input range: [-infinity , +infinity).
		 Search range: [-1e4, 1E4]

      STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		 1 if answer appears to be lower than lowest
		   search bound
		 2 if answer appears to be higher than greatest
		   search bound
		 3 if P + Q != 1

      BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.

				 Method

      Upper tail    of	the  cumulative	 noncentral t is calculated usin
      formulae	from page 532  of Johnson, Kotz,  Balakrishnan, Coninuou
      Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)

      Computation of other parameters involve a seach for a value that
      produces	the desired  value  of P.   The search relies  on  the
      monotinicity of P with the other parameter.

**********************************************************************/
#define tent4 1e4
#define tol 1e-8
#define atol 1e-50
#define Azero 1e-100
#define one ( 1 - 1e-16 )
#define inf 1e100

    static double K3 = 0.5;
    static double K4 = 5;
    static double ccum,cum,fx;
    static logical qhi,qleft;
    static double T1,T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14;
    /*
      ..
      .. Executable Statements ..
    */

    DBGprt5("cdftnc(%d, p=%g, t=%g, df=%g", *which, *p, *t, *df);
    DBGprt2(", ncp= %g)\n", *pnonc);

    which_check_1(4);
    I001_check(1,*p,-2);
    I001_check(1,*q,-3);

    if(*which != 3 &&  *df  <=	0) { *bound = 0; *status = -5; return; }

/*---- regular :  which in {1,2,3,4},	p & df in range : */

    switch(*which) {
    case 1:
	cumtnc(t,df,pnonc,p,q);
	*status = 0;
	return;

    case 2:
	*t = 5;
	T1 = -inf;
	T2 = inf;
	T5 = atol;
	T6 = tol;
	dstinv(&T1,&T2,&K3,&K3,&K4,&T5,&T6);
	*status = 0;
	dinvr(status,t,&fx,&qleft,&qhi);

	while(*status == 1) {/* inversion iterations */
	    cumtnc(t,df,pnonc,&cum,&ccum);
	    fx = cum - *p;
	    dinvr(status,t,&fx,&qleft,&qhi);
	}

	if(*status == -1){
	    if(qleft) {
		*status = 1; *bound = -inf;
	    } else {
		*status = 2; *bound = inf;
	    }
	}
	return;

    case 3:
	*df = 5;
	T7 = Azero;
	T8 = tent4;
	T9 = atol;
	T10 = tol;
	dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
	*status = 0;
	dinvr(status,df,&fx,&qleft,&qhi);

	while(*status == 1){
	    cumtnc(t,df,pnonc,&cum,&ccum);
	    fx = cum - *p;
	    dinvr(status,df,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = Azero;
	    } else {
		*status = 2; *bound = inf;
	    }
	}
	return;

    case 4:
	*pnonc = 5;
	T11 = -tent4;
	T12 = tent4;
	T13 = atol;
	T14 = tol;
	dstinv(&T11,&T12,&K3,&K3,&K4,&T13,&T14);
	*status = 0;
	dinvr(status,pnonc,&fx,&qleft,&qhi);

	while(*status == 1) {
	    cumtnc(t,df,pnonc,&cum,&ccum);
	    fx = cum - *p;
	    dinvr(status,pnonc,&fx,&qleft,&qhi);
	}
	if(*status == -1) {
	    if(qleft) {
		*status = 1; *bound = 0;
	    } else {
		*status = 2; *bound = tent4;
	    }
	}
	return;
    }
#undef tent4
#undef tol
#undef atol
#undef Azero
#undef one
#undef inf
} /* cdftnc() */

void cumbet(double *x,double *y,double *a,double *b,double *cum,double *ccum)
{
/*
**********************************************************************

     void cumbet(double *x,double *y,double *a,double *b,double *cum,
	    double *ccum)

	  Double precision cUMulative incomplete BETa distribution


			      Function


     Calculates the cdf to X of the incomplete beta distribution
     with parameters a and b.  This is the integral from 0 to x
     of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)


			      Arguments


     X --> Upper limit of integration.
					X is DOUBLE PRECISION

     Y --> 1 - X.
					Y is DOUBLE PRECISION

     A --> First parameter of the beta distribution.
					A is DOUBLE PRECISION

     B --> Second parameter of the beta distribution.
					B is DOUBLE PRECISION

     CUM <-- Cumulative incomplete beta distribution.
					CUM is DOUBLE PRECISION

     CCUM <-- Compliment of Cumulative incomplete beta distribution.
					CCUM is DOUBLE PRECISION


			      Method


     Calls the routine BRATIO.

				   References

     Didonato, Armido R. and Morris, Alfred H. Jr. (1992) Algorithim
     708 Significant Digit Computation of the Incomplete Beta Function
     Ratios. ACM ToMS, Vol.18, No. 3, Sept. 1992, 360-373.

**********************************************************************
*/
    int ierr;

    if(*x <= 0) {
	*cum = 0; *ccum = 1; return;
    }
    if(*y <= 0) {
	*cum = 1; *ccum = 0; return;
    }
    bratio(a,b,x,y,cum,ccum,&ierr);
    return;
} /* cumbet() */

void cumbin(double *s,double *xn,double *pr,double *ompr,
	    double *cum,double *ccum)
{
/*
**********************************************************************

     void cumbin(double *s,double *xn,double *pr,double *ompr,
	    double *cum,double *ccum)

		    CUmulative BINomial distribution


			      Function


     Returns the probability   of 0  to	 S  successes in  XN   binomial
     trials, each of which has a probability of success, PBIN.


			      Arguments


     S --> The upper limit of cumulation of the binomial distribution.
						  S is DOUBLE PRECISION

     XN --> The number of binomial trials.
						  XN is DOUBLE PRECISIO

     PBIN --> The probability of success in each binomial trial.
						  PBIN is DOUBLE PRECIS

     OMPR --> 1 - PBIN
						  OMPR is DOUBLE PRECIS

     CUM <-- Cumulative binomial distribution.
						  CUM is DOUBLE PRECISI

     CCUM <-- Compliment of Cumulative binomial distribution.
						  CCUM is DOUBLE PRECIS


			      Method


     Formula  26.5.24	 of   Abramowitz  and	 Stegun,  Handbook   of
     Mathematical   Functions (1966) is	  used	to reduce the  binomial
     distribution  to  the  cumulative	  beta distribution.

**********************************************************************
*/
    double T1,T2;

    if(*s < *xn) {
	T1 = *s+1;
	T2 = *xn-*s;
	cumbet(pr,ompr,&T1,&T2,ccum,cum);
    } else {
	*cum = 1;
	*ccum = 0;
    }
    return;
} /* cumbin() */

void cumchi(double *x,double *df,double *cum,double *ccum)
{
/*
**********************************************************************

     void cumchi(double *x,double *df,double *cum,double *ccum)
	     CUMulative of the CHi-square distribution


			      Function


     Calculates the cumulative chi-square distribution.


			      Arguments


     X	     --> Upper limit of integration of the
		 chi-square distribution.
						 X is DOUBLE PRECISION

     DF	     --> Degrees of freedom of the
		 chi-square distribution.
						 DF is DOUBLE PRECISION

     CUM <-- Cumulative chi-square distribution.
						 CUM is DOUBLE PRECISIO

     CCUM <-- Compliment of Cumulative chi-square distribution.
						 CCUM is DOUBLE PRECISI


			      Method


     Calls incomplete gamma function (CUMGAM)

**********************************************************************
*/
    double a,xx;

    a = *df*0.5;
    xx = *x*0.5;
    cumgam(&xx,&a,cum,ccum);
    return;
} /* cumchi() */

void cumchn(double *x,double *df,double *pnonc, double *cum, double *ccum)
{
/**********************************************************************

     void cumchn(double *x,double *df,double *pnonc,double *cum,
		 double *ccum)

	     CUMulative of the Non-central CHi-square distribution

			       Function

     Calculates	    the	      cumulative      non-central    chi-square
     distribution, i.e.,  the probability   that  a   random   variable
     which    follows  the  non-central chi-square  distribution,  with
     non-centrality  parameter	  PNONC	 and   continuous  degrees   of
     freedom DF, is less than or equal to X.

			      Arguments

     X	     --> Upper limit of integration of the non-central
		 chi-square distribution.

     DF	     --> Degrees of freedom of the non-central
		 chi-square distribution.

     PNONC   --> Non-centrality parameter of the non-central
		 chi-square distribution.

     CUM <-- Cumulative non-central chi-square distribution.

     CCUM <-- Compliment of Cumulative non-central chi-square distribution.


				Method

     Uses  formula  26.4.25   of  Abramowitz  and  Stegun, Handbook  of
     Mathematical    Functions,	 US   NBS   (1966)    to calculate  the
     non-central chi-square.

				Variables

     EPS     --- Convergence criterion.	 The sum stops when a
		 term is less than EPS*SUM
		 _or_  SUM < 1e-20  [i.e. have relative & absolute conv.check]

**********************************************************************/
#define dg(i) (*df + 2 * (double)(i))

static const double eps = 1e-5;

    double adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt,lfact,pcent,pterm,sum,
	sumadj,term,wt,xnonc;
    int i,icent;
    double T1,T2,T3;

    if(*x <= 0) {
	*cum = 0; *ccum = 1; return;
    }

    if(*pnonc <= 1e-10) {
	/*
	  When non-centrality parameter is (essentially) zero,
	  use cumulative chi-square distribution
	*/
	cumchi(x,df,cum,ccum);
	return;
    }

    xnonc = *pnonc / 2;
/*
***********************************************************************
     The following code calcualtes the weight, chi-square, and
     adjustment term for the central term in the infinite series.
     The central term is the one in which the poisson weight is
     greatest.	The adjustment term is the amount that must
     be subtracted from the chi-square to move up two degrees
     of freedom.
***********************************************************************
*/
    icent = fifidint(xnonc);
    if(icent == 0) icent = 1;
    chid2 = *x / 2;
/*
     Calculate central weight term
*/
    T1 = (double)(icent + 1);
    lfact = alngam(&T1);
    lcntwt = -xnonc + (double)icent * log(xnonc) - lfact;
    centwt = exp(lcntwt);
/*
     Calculate central chi-square
*/
    T2 = dg(icent);
    cumchi(x,&T2,&pcent,ccum);
/*
     Calculate central adjustment term
*/
    dfd2 = dg(icent) / 2;
    T3 = 1 + dfd2;
    lfact = alngam(&T3);
    lcntaj = dfd2 * log(chid2) - chid2 - lfact;
    centaj = exp(lcntaj);
    sum = centwt * pcent;
/*
***********************************************************************
     Sum backwards from the central term towards zero.
     Quit whenever either
     (1) the zero term is reached, or
     (2) the term gets small relative to the sum
***********************************************************************
*/
    sumadj = 0;
    adj = centaj;
    wt = centwt;
    i = icent;

    do {
	dfd2 = dg(i) / 2;
	/* Adjust chi-square for two fewer degrees of freedom.
	   The adjusted value ends up in PTERM.
	*/
	adj = adj * dfd2 / chid2;
	sumadj += adj;
	pterm = pcent + sumadj;
	/* Adjust poisson weight for J decreased by one */
	wt *= ((double)i / xnonc);
	term = wt * pterm;
	sum += term;
	i -= 1;
    } while(!(qsmall(term, sum, 1e-20, eps) || i == 0));
/*
***********************************************************************
     Now sum forward from the central term towards infinity.
     Quit when either
     (1) the term gets small relative to the sum, or
***********************************************************************
*/
    sumadj = adj = centaj;
    wt = centwt;
    i = icent;
    do {
	/* Update weights for next higher J */
	wt *= (xnonc / (double)(i + 1));
	/* Calculate PTERM and add term to sum */
	pterm = pcent - sumadj;
	term = wt * pterm;
	sum += term;
	/* Update adjustment term for DF for next iteration */
	i += 1;
	dfd2 = dg(i) / 2;
	adj = adj * chid2 / dfd2;
	sumadj += adj;
    } while(!qsmall(term, sum, 1e-20, eps));
    *cum = sum;
    *ccum = 0.5 + (0.5 - *cum);
    return;
#undef dg
} /* cumchn() */

void cumf(double *f,double *dfn,double *dfd,double *cum,double *ccum)
{
/*
**********************************************************************

     void cumf(double *f,double *dfn,double *dfd,double *cum,double *ccum)
		    CUMulative F distribution


			      Function


     Computes  the  integral from  0  to  F of	the f-density  with DFN
     and DFD degrees of freedom.


			      Arguments


     F --> Upper limit of integration of the f-density.
						  F is DOUBLE PRECISION

     DFN --> Degrees of freedom of the numerator sum of squares.
						  DFN is DOUBLE PRECISI

     DFD --> Degrees of freedom of the denominator sum of squares.
						  DFD is DOUBLE PRECISI

     CUM <-- Cumulative f distribution.
						  CUM is DOUBLE PRECISI

     CCUM <-- Compliment of Cumulative f distribution.
						  CCUM is DOUBLE PRECIS


			      Method


     Formula  26.5.28 of  Abramowitz and   Stegun   is	used to	 reduce
     the cumulative F to a cumulative beta distribution.


			      Note


     If F is less than or equal to 0, 0 is returned.

**********************************************************************
*/
#define half 0.5
#define one 1
    double dsum,prod,xx,yy;
    int ierr;
    double T1,T2;

    if(*f <= 0) {
	*cum = 0; *ccum = 1; return;
    }

    prod = *dfn**f;
/*
     XX is such that the incomplete beta with parameters
     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
     YY is 1 - XX
     Calculate the smaller of XX and YY accurately
*/
    dsum = *dfd+prod;
    xx = *dfd/dsum;
    if(xx > half) {
	yy = prod/dsum;
	xx = one-yy;
    }
    else
	yy = one-xx;
    T1 = *dfd*half;
    T2 = *dfn*half;
    bratio(&T1,&T2,&xx,&yy,ccum,cum,&ierr);
    return;
#undef half
#undef one
} /* cumf */

void cumfnc(double *f,double *dfn,double *dfd,double *pnonc,
	    double *cum,double *ccum)
{
/*
**********************************************************************

	       F -NON- -C-ENTRAL F DISTRIBUTION



			      Function


     COMPUTES NONCENTRAL F DISTRIBUTION WITH DFN AND DFD
     DEGREES OF FREEDOM AND NONCENTRALITY PARAMETER PNONC


			      Arguments


     X --> UPPER LIMIT OF INTEGRATION OF NONCENTRAL F IN EQUATION

     DFN --> DEGREES OF FREEDOM OF NUMERATOR

     DFD -->  DEGREES OF FREEDOM OF DENOMINATOR

     PNONC --> NONCENTRALITY PARAMETER.

     CUM <-- CUMULATIVE NONCENTRAL F DISTRIBUTION

     CCUM <-- COMPLIMENT OF CUMMULATIVE


			      Method


     USES FORMULA 26.6.20 OF REFERENCE FOR INFINITE SERIES.
     SERIES IS CALCULATED BACKWARD AND FORWARD FROM J = LAMBDA/2
     (THIS IS THE TERM WITH THE LARGEST POISSON WEIGHT) UNTIL
     THE CONVERGENCE CRITERION IS MET.

     FOR SPEED, THE INCOMPLETE BETA FUNCTIONS ARE EVALUATED
     BY FORMULA 26.5.16.


	       REFERENCE


     HANDBOOD OF MATHEMATICAL FUNCTIONS
     EDITED BY MILTON ABRAMOWITZ AND IRENE A. STEGUN
     NATIONAL BUREAU OF STANDARDS APPLIED MATEMATICS SERIES - 55
     MARCH 1965
     P 947, EQUATIONS 26.6.17, 26.6.18
*/

/*
			      Note

     The sum continues until a succeeding term is less than eps
     times the sum  _or_ the sum is less than 1e-20).
     [i.e. have relative & absolute conv.check]

     eps is set to 1e-4 in a data statement which can be changed.
*/
static const double eps = 1e-4;

#define half 0.5
#define one 1

    double dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,
	upterm,xmult,xnonc;
    int i,icent,ierr;
    double T1,T2,T3,T4,T5,T6;

    if(*f <= 0) {
	*cum = 0; *ccum = 1; return;
    }

    if(*pnonc < 1e-10) {/* The non-centrality parameter is (essentially) zero.
			 */
	cumf(f,dfn,dfd,cum,ccum);
	return;
    }

    xnonc = *pnonc/2;
/*
     Calculate the central term of the poisson weighting factor.
*/
    icent = (long)(xnonc);
    if(icent == 0) icent = 1;
/*
     Compute central weight term
*/
    T1 = (double)(icent+1);
    centwt = exp(-xnonc+(double)icent*log(xnonc)-alngam(&T1));
/*
     Compute central incomplete beta term
     Assure that minimum of arg to beta and 1 - arg is computed
	  accurately.
*/
    prod = *dfn**f;
    dsum = *dfd+prod;
    yy = *dfd/dsum;
    if(yy > half) {
	xx = prod/dsum;
	yy = one-xx;
    }
    else
	xx = one-yy;
    T2 = *dfn*half+(double)icent;
    T3 = *dfd*half;
    bratio(&T2,&T3,&xx,&yy,&betdn,&dummy,&ierr);
    adn = *dfn/2+(double)icent;
    aup = adn;
    b = *dfd/2;
    betup = betdn;
    sum = centwt*betdn;
/*
     Now sum terms backward from icent until convergence or all done
*/
    xmult = centwt;
    i = icent;
    T4 = adn+b;
    T5 = adn+1;
    dnterm = exp(alngam(&T4)-alngam(&T5)-alngam(&b)+adn*log(xx)+b*log(yy));

    while(!qsmall(xmult*betdn, sum, 1e-20, eps) && i > 0) {
	xmult *= ((double)i/xnonc);
	i -= 1;
	adn -= 1;
	dnterm = (adn+1)/((adn+b)*xx)*dnterm;
	betdn += dnterm;
	sum += (xmult*betdn);
    }

    i = icent+1;
/*
     Now sum forwards until convergence
*/
    xmult = centwt;
    if(aup-1+b == 0)
	upterm = exp( -alngam(&aup) - alngam(&b)+ (aup-1)*log(xx)+ b*log(yy));
    else  {
	T6 = aup-1+b;
	upterm = exp(alngam(&T6)
		      -alngam(&aup) - alngam(&b)+ (aup-1)*log(xx)+ b*log(yy));
    }

    do {
	xmult *= (xnonc/(double)i);
	i += 1;
	aup += 1;
	upterm = (aup+b-2)*xx/(aup-1)*upterm;
	betup -= upterm;
	sum += (xmult*betup);
    } while(!qsmall(xmult*betup, sum, 1e-20, eps));

    *cum = sum;
    *ccum = 0.5+(0.5-*cum);
    return;
#undef half
#undef one
} /* cumfnc() */

void cumgam(double *x,double *a,double *cum,double *ccum)
{
/*
**********************************************************************

     void cumgam(double *x,double *a,double *cum,double *ccum)
	   Double precision cUMulative incomplete GAMma distribution


			      Function


     Computes	the  cumulative	       of    the     incomplete	  gamma
     distribution, i.e., the integral from 0 to X of
	  (1/GAM(A))*EXP(-T)*T**(A-1) DT
     where GAM(A) is the complete gamma function of A, i.e.,
	  GAM(A) = integral from 0 to infinity of
		    EXP(-T)*T**(A-1) DT


			      Arguments


     X --> The upper limit of integration of the incomplete gamma.
						X is DOUBLE PRECISION

     A --> The shape parameter of the incomplete gamma.
						A is DOUBLE PRECISION

     CUM <-- Cumulative incomplete gamma distribution.
					CUM is DOUBLE PRECISION

     CCUM <-- Compliment of Cumulative incomplete gamma distribution.
						CCUM is DOUBLE PRECISIO


			      Method


     Calls the routine GRATIO.

**********************************************************************
*/
    static int K0 = 0;

    if(*x <= 0) {
	*cum = 0; *ccum = 1; return;
    }
    gratio(a,x,cum,ccum,&K0);
    return;
} /* cumgam() */

void cumnbn(double *s,double *xn,double *pr,double *ompr,
	    double *cum,double *ccum)
{
/*
**********************************************************************

     void cumnbn(double *s,double *xn,double *pr,double *ompr,
	    double *cum,double *ccum)

		    CUmulative Negative BINomial distribution


			      Function


     Returns the probability that it there will be S or fewer failures
     before there are XN successes, with each binomial trial having
     a probability of success PR.

     Prob(# failures = S | XN successes, PR)  =
			( XN + S - 1 )
			(	     ) * PR^XN * (1-PR)^S
			(      S     )


			      Arguments


     S --> The number of failures
						  S is DOUBLE PRECISION

     XN --> The number of successes
						  XN is DOUBLE PRECISIO

     PR --> The probability of success in each binomial trial.
						  PR is DOUBLE PRECISIO

     OMPR --> 1 - PR
						  OMPR is DOUBLE PRECIS

     CUM <-- Cumulative negative binomial distribution.
						  CUM is DOUBLE PRECISI

     CCUM <-- Compliment of Cumulative negative binomial distribution.
						  CCUM is DOUBLE PRECIS


			      Method


     Formula  26.5.26	 of   Abramowitz  and	 Stegun,  Handbook   of
     Mathematical   Functions (1966) is	  used	to reduce the  negative
     binomial distribution to the cumulative beta distribution.

**********************************************************************
*/
    double T1;

    T1 = *s+1.e0;
    cumbet(pr,ompr,xn,&T1,cum,ccum);
    return;
} /* cumnbn() */

void cumnor(double *arg,double *result,double *ccum)
{
/*
**********************************************************************

     void cumnor(double *arg,double *result,double *ccum)


			      Function


     Computes the cumulative  of    the	 normal	  distribution,	  i.e.,
     the integral from -infinity to x of
	  (1/sqrt(2*pi)) exp(-u*u/2) du

     X --> Upper limit of integration.
					X is DOUBLE PRECISION

     RESULT <-- Cumulative normal distribution.
					RESULT is DOUBLE PRECISION

     CCUM <-- Compliment of Cumulative normal distribution.
					CCUM is DOUBLE PRECISION

     Renaming of function ANORM from:

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     with slight modifications to return ccum and to deal with
     machine constants.

**********************************************************************
  Original Comments:
------------------------------------------------------------------

 This function evaluates the normal distribution function:

			      / x
		     1	     |	     -t*t/2
	  P(x) = ----------- |	    e	    dt
		 sqrt(2 pi)  |
			     /-oo

   The main computation evaluates near-minimax approximations
   derived from those in "Rational Chebyshev approximations for
   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
   This transportable program uses rational functions that
   theoretically approximate the normal distribution function to
   at least 18 significant decimal digits.  The accuracy achieved
   depends on the arithmetic system, the compiler, the intrinsic
   functions, and proper selection of the machine-dependent
   constants.

*******************************************************************
*******************************************************************

 Explanation of machine-dependent constants.

   MIN	 = smallest machine representable number.

   EPS	 = argument below which anorm(x) may be represented by
	   0.5	and above which	 x*x  will not underflow.
	   A conservative value is the largest machine number X
	   such that   1 + X = 1   to machine precision.
*******************************************************************
*******************************************************************

 Error returns

  The program returns  ANORM = 0     for  ARG <= XLOW.


 Intrinsic functions required are:

     ABS, AINT, EXP


  Author: W. J. Cody
	  Mathematics and Computer Science Division
	  Argonne National Laboratory
	  Argonne, IL 60439

  Latest modification: March 15, 1992

------------------------------------------------------------------
*/
static double a[5] = {
    2.2352520354606839287,1.6102823106855587881e02,1.0676894854603709582e03,
    1.8154981253343561249e04,6.5682337918207449113e-2
};
static double b[4] = {
    4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
    4.5507789335026729956e04
};
static double c[9] = {
    3.9894151208813466764e-1,8.8831497943883759412,9.3506656132177855979e01,
    5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
    1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
};
static double d[8] = {
    2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
    6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
    3.8912003286093271411e04,1.9685429676859990727e04
};
static double half = 0.5;
static double p[6] = {
    2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
    1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
};
static double one = 1;
static double q[5] = {
    1.28426009614491121,4.68238212480865118e-1,6.59881378689285515e-2,
    3.78239633202758244e-3,7.29751555083966205e-5
};
static double sixten = 1.60;
static double sqrpi = 3.9894228040143267794e-1;
static double thrsh = 0.66291;
static double root32 = 5.656854248;
static double zero = 0;

    int i;
    double del,eps,temp,x,xden,xnum,y,xsq,min;
/*
------------------------------------------------------------------
  Machine dependent constants
------------------------------------------------------------------
*/
    eps = spmpar(1)*0.5;
    min = spmpar(2);
    x = *arg;
    y = fabs(x);
    if(y <= thrsh) {
	/*
	  ------------------------------------------------------------------
	  Evaluate  anorm  for	|X| <= 0.66291
	  ------------------------------------------------------------------
	*/
	xsq = zero;
	if(y > eps) xsq = x*x;
	xnum = a[4]*xsq;
	xden = xsq;
	for(i=0; i<3; i++) {
	    xnum = (xnum+a[i])*xsq;
	    xden = (xden+b[i])*xsq;
	}
	*result = x*(xnum+a[3])/(xden+b[3]);
	temp = *result;
	*result = half+temp;
	*ccum = half-temp;
    }
    else if(y <= root32) {
	/*
	  ------------------------------------------------------------------
	  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
	  ------------------------------------------------------------------
	*/
	xnum = c[8]*y;
	xden = y;
	for(i=0; i<7; i++) {
	    xnum = (xnum+c[i])*y;
	    xden = (xden+d[i])*y;
	}
	*result = (xnum+c[7])/(xden+d[7]);
	xsq = fifdint(y*sixten)/sixten;
	del = (y-xsq)*(y+xsq);
	*result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
	*ccum = one-*result;
	if(x > zero) {
	    temp = *result;
	    *result = *ccum;
	    *ccum = temp;
	}
    }
    else  {
	/*
	  ------------------------------------------------------------------
	  Evaluate  anorm  for |X| > sqrt(32)
	  ------------------------------------------------------------------
	*/
	*result = zero;
	xsq = one/(x*x);
	xnum = p[5]*xsq;
	xden = xsq;
	for(i=0; i<4; i++) {
	    xnum = (xnum+p[i])*xsq;
	    xden = (xden+q[i])*xsq;
	}
	*result = xsq*(xnum+p[4])/(xden+q[4]);
	*result = (sqrpi-*result)/y;
	xsq = fifdint(x*sixten)/sixten;
	del = (x-xsq)*(x+xsq);
	*result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
	*ccum = one-*result;
	if(x > zero) {
	    temp = *result;
	    *result = *ccum;
	    *ccum = temp;
	}
    }
    if(*result < min) *result = 0;
/*
------------------------------------------------------------------
  Fix up for negative argument, erf, etc.
------------------------------------------------------------------
*/
    if(*ccum < min) *ccum = 0;
} /* cumnor() */

void cumpoi(double *s,double *xlam,double *cum,double *ccum)
{
/*
**********************************************************************

     void cumpoi(double *s,double *xlam,double *cum,double *ccum)
		    CUMulative POIsson distribution


			      Function


     Returns the  probability  of  S   or  fewer events in  a	Poisson
     distribution with mean XLAM.


			      Arguments


     S --> Upper limit of cumulation of the Poisson.
						  S is DOUBLE PRECISION

     XLAM --> Mean of the Poisson distribution.
						  XLAM is DOUBLE PRECIS

     CUM <-- Cumulative poisson distribution.
					CUM is DOUBLE PRECISION

     CCUM <-- Compliment of Cumulative poisson distribution.
						  CCUM is DOUBLE PRECIS


			      Method


     Uses formula  26.4.21   of	  Abramowitz and  Stegun,  Handbook  of
     Mathematical   Functions  to reduce   the	 cumulative Poisson  to
     the cumulative chi-square distribution.

**********************************************************************
*/
    double chi,df;

    df = 2*(*s+1);
    chi = 2**xlam;
    cumchi(&chi,&df,ccum,cum);
    return;
} /* cumpoi() */

void cumt(double *t,double *df,double *cum,double *ccum)
{
/*
**********************************************************************

     void cumt(double *t,double *df,double *cum,double *ccum)
		    CUMulative T-distribution


			      Function


     Computes the integral from -infinity to T of the t-density.


			      Arguments


     T --> Upper limit of integration of the t-density.
						  T is DOUBLE PRECISION

     DF --> Degrees of freedom of the t-distribution.
						  DF is DOUBLE PRECISIO

     CUM <-- Cumulative t-distribution.
						  CCUM is DOUBLE PRECIS

     CCUM <-- Compliment of Cumulative t-distribution.
						  CCUM is DOUBLE PRECIS


			      Method


     Formula 26.5.27   of     Abramowitz  and	Stegun,	   Handbook  of
     Mathematical Functions  is	  used	 to  reduce the	 t-distribution
     to an incomplete beta.

**********************************************************************
*/
    static double Half = 0.5;
    static double xx,a,oma,tt,yy,dfptt,T1;

    tt = *t**t;
    dfptt = *df+tt;
    xx = *df/dfptt;/*  nu/(nu + t^2) */
    yy = tt/dfptt; /* t^2/(nu + t^2) */
    T1 = 0.5**df;
    cumbet(&xx,&yy,&T1,&Half,&a,&oma);
    if(*t <= 0) {
	*cum = 0.5*a;
	*ccum = oma+*cum;
    }
    else {
	*ccum = 0.5*a;
	*cum = oma+*ccum;
    }
    return;
} /* cumt() */

void cumtnc(double *t,double *df,double *pnonc,
	    double *cum, double *ccum)
{
/**********************************************************************

     void cumtnc(double *t,double *df,double *pnonc,double *cum,
		 double *ccum)

		  CUMulative Non-Central T-distribution


			       Function


      Computes the integral from -infinity to T of the non-central
      t-density.


			       Arguments


      T --> Upper limit of integration of the non-central t-density.

      DF --> Degrees of freedom of the non-central t-distribution.

      PNONC --> Non-centrality parameter of the non-central t distibutio

      CUM <-- Cumulative t-distribution.

      CCUM <-- Compliment of Cumulative t-distribution.


			       Method

      Upper tail    of	the  cumulative	 noncentral t	using
      formulae from page 532  of Johnson, Kotz,	 Balakrishnan, Coninuous
      Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)

      This implementation starts the calculation at i = lambda,
      which is near the largest Di.  It then sums forward and backward.
**********************************************************************/
#define one 1
#define zero 0
#define half 0.5
#define two 2
#define onep5 1.5
#define conv 1e-7/* convergence tolerance */
#define tiny 1e-10
    static double alghdf,b,bb,bbcent,bcent,cent,d,dcent,dpnonc,dum1,dum2,
	e,ecent,halfdf,lambda,lnomx,lnx,omx,pnonc2,s,scent,ss,sscent,
	t2,term,tt,twoi,x,xi, xlnd,xlne;
    static int ierr;
    static logical qrevs;
    static double T1,T2,T3,T4,T5,T6,T7,T8,T9,T10;

    if(fabs(*pnonc) <= tiny) { /* pnonc essentially zero: use (central) t */
	cumt(t,df,cum,ccum);
	return;
    }
    qrevs = *t < zero;
    if(qrevs) {
	tt = -*t;
	dpnonc = -*pnonc;
    }
    else  {
	tt = *t;
	dpnonc = *pnonc;
    }
    t2 = tt * tt;
    if(fabs(tt) <= tiny) {/* t ~ 0 ==> Phi( -delta ) */
	T1 = -*pnonc;
	cumnor(&T1,cum,ccum);
	return;
    }
    pnonc2 = dpnonc * dpnonc;
    lambda = half * pnonc2;
    x = *df / (*df + t2);
    omx = one - x;
    lnx = log(x);
    lnomx = log(omx);
    halfdf = half * *df;
    alghdf = gamln(&halfdf);
/*
******************** Case i = lambda
*/
    cent = fifidint(lambda);
    if(cent < one) cent = one;
/*
  Compute d=T(2i) in log space and offset by exp(-lambda)
*/
    T2 = cent + one;
    xlnd = cent * log(lambda) - gamln(&T2) - lambda;
    dcent = exp(xlnd);
/*
  Compute e=t(2i+1) in log space offset by exp(-lambda)
*/
    T3 = cent + onep5;
    xlne = (cent + half) * log(lambda) - gamln(&T3) - lambda;
    ecent = exp(xlne);
    if(dpnonc < zero) ecent = -ecent;
/*
  Compute bcent=B(2*cent)
*/
    T4 = cent + half;
    bratio(&halfdf,&T4,&x,&omx,&bcent,&dum1,&ierr);
/*
  compute bbcent=B(2*cent+1)
*/
    T5 = cent + one;
    bratio(&halfdf,&T5,&x,&omx,&bbcent,&dum2,&ierr);

    if(bcent + bbcent < tiny) {
	/* bcent and bbcent are essentially zero
	   Thus t is effectively infinite
	*/
	if(qrevs) {
	    *cum = zero;
	    *ccum = one;
	}
	else  {
	    *cum = one;
	    *ccum = zero;
	}
	return;
    }

    if(dum1 + dum2 < tiny) {
	/* bcent and bbcent are essentially one
	   Thus t is effectively zero
	*/
	T6 = -*pnonc;
	cumnor(&T6,cum,ccum);
	return;
    }
/*
  First term in ccum is D*B + E*BB
*/
    *ccum = dcent * bcent + ecent * bbcent;
/*
  compute s(cent) = B(2*(cent+1)) - B(2*cent))
*/
    T7 = halfdf + cent + half;
    T8 = cent + onep5;
    scent = gamln(&T7) - gamln(&T8) - alghdf + halfdf * lnx + (cent + half) *
	lnomx;
    scent = exp(scent);
/*
  compute ss(cent) = B(2*cent+3) - B(2*cent+1)
*/
    T9 = halfdf + cent + one;
    T10 = cent + two;
    sscent = gamln(&T9) - gamln(&T10) - alghdf + halfdf * lnx + (cent + one) *
	lnomx;
    sscent = exp(sscent);
/*
******************** Sum Forward
*/
    xi = cent + one;
    twoi = two * xi;
    d = dcent;
    e = ecent;
    b = bcent;
    bb = bbcent;
    s  = scent;
    ss = sscent;

    do {
	b += s;
	bb += ss;
	d = lambda / xi * d;
	e = lambda / (xi + half) * e;
	term = d * b + e * bb;
	*ccum += term;
	s  *= omx * (*df + twoi - one) / (twoi + one);
	ss *= omx * (*df + twoi) / (twoi + two);
	xi += one;
	twoi = two * xi;
    } while(fabs(term) > conv * *ccum);
/*
******************** Sum Backward
*/
    xi = cent;
    twoi = two * xi;
    d = dcent;
    e = ecent;
    b = bcent;
    bb = bbcent;
    s  = scent	* (one + twoi) / ((*df + twoi - one) * omx);
    ss = sscent * (two + twoi) / ((*df + twoi) * omx);

    do {
	b -= s;
	bb -= ss;
	d *= (xi / lambda);
	e *= ((xi + half) / lambda);
	term = d * b + e * bb;
	*ccum += term;
	xi -= one;
	if(xi < half) break;
	twoi = two * xi;
	s = s * (one + twoi) / ((*df + twoi - one) * omx);
	ss = ss * (two + twoi) / ((*df + twoi) * omx);
    } while(fabs(term) > conv * *ccum);

    if(qrevs) {
	*cum = half * *ccum;
	*ccum = one - *cum;
    }
    else  {
	*ccum = half * *ccum;
	*cum = one - *ccum;
    }
/*
  Due to roundoff error the answer may not lie between zero and one
  Force it to do so
*/
    *cum  = fmax2(fmin2(*cum, one),zero);
    *ccum = fmax2(fmin2(*ccum,one),zero);
    return;
#undef one
#undef zero
#undef half
#undef two
#undef onep5
#undef conv
#undef tiny
} /* cumtnc() */

#undef which_check_1
#undef I01_check
#undef I001_check
#undef PpQ_check
#undef qsmall


double devlpl(double a[], int n, double x)
{
/*
**********************************************************************

     double devlpl(double a[],int *n,double *x)

		      EVALuate a PoLynomial at X

			      Function

     returns
	  A(1) + A(2)*X + ... + A(N)*X**(N-1)


			      Arguments


     A --> Array of coefficients of the polynomial.
					A is DOUBLE PRECISION(N)

     N --> Length of A, also degree of polynomial - 1.
					N is INTEGER

     X --> Point at which the polynomial is to be evaluated.
					X is DOUBLE PRECISION

**********************************************************************
*/
    double term;
    int i;

    term = a[n-1];
    for(i= n-1-1; i>=0; i--)
	term = a[i] + term*x;
    return term;
}


double dinvnr(double *p,double *q)
{
/*
**********************************************************************

     double dinvnr(double *p,double *q)
     Double precision NoRmal distribution INVerse  -- aka "qnorm()"


			      Function


     Returns X	such that CUMNOR(X)  =	 P,  i.e., the	integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


			      Arguments


     P --> The probability whose normal deviate is sought.
		    P is DOUBLE PRECISION

     Q --> 1-P
		    P is DOUBLE PRECISION


			      Method


     The  rational   function	on  page 95    of Kennedy  and	Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
     value for the Newton method of finding roots.


			      Note


     If P or Q < machine EPS returns +/- DINVNR(EPS)

**********************************************************************
*/
#define maxit 100
#define eps 1e-13
#define nhalf -0.5
#define dennor(x) (M_1_SQRT_2PI*exp(nhalf*(x)*(x)))
static double dinvnr,strtx,xcur,cum,ccum,pp,dx;
static int i;
static logical qporq;
/*
     ..
     .. Executable Statements ..
*/
/*
     FIND MINIMUM OF P AND Q
*/
    qporq = *p <= *q;
    if(!qporq) goto S10;
    pp = *p;
    goto S20;
S10:
    pp = *q;
S20:
/*
     INITIALIZATION STEP
*/
    strtx = stvaln(&pp);
    xcur = strtx;
/*
     NEWTON INTERATIONS
*/
    for(i=1; i<=maxit; i++) {
	cumnor(&xcur,&cum,&ccum);
	dx = (cum-pp)/dennor(xcur);
	xcur -= dx;
	if(fabs(dx/xcur) < eps) goto S40;
    }
    dinvnr = strtx;
/*
	Newton has failed :
*/
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
S40:
/*
	Newton has succeeded :
*/
    dinvnr = xcur;
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
#undef maxit
#undef eps
#undef nhalf
#undef dennor
} /* dinvnr() */

/* DEFINE DINVR */
static void E0000(int IENTRY,
		  /* those for	dinvr() : */
		  int *status,double *x,double *fx,
		  logical *qleft,logical *qhi,
		  /* these for	dstinv() : */
		  double *zabsst,
		  double *zabsto,double *zbig,double *zrelst,
		  double *zrelto,double *zsmall,double *zstpmu)
{
#define is_monotone(x,y,z) (int)((x) <= (y) && (y) <= (z))

    /* Globals initialized from	 dstinv() : */
    static double absstp,abstol, relstp,reltol, big, small,stpmul;

    static double fbig,fsmall,step,xhi,xlb,xlo,xsave,xub,yy;
    static int i99999;
    static logical qbdd,qcond,qdum1,qdum2,qincr,qlim,qok,qup;

    if(IENTRY == 1) { /*--- dstinv() : Initialize --- */
	small = *zsmall;
	big = *zbig;
	absstp = *zabsst;
	relstp = *zrelst;
	stpmul = *zstpmu;
	abstol = *zabsto;
	reltol = *zrelto;
	return;
    }

    /* Else --- DINVR --- */

    if(*status <= 0) {

	if(!is_monotone(small,*x,big))
	    ftnstop("(SMALL, X, BIG) not monotone in INVR");
	xsave = *x;
	/*
	  See that SMALL and BIG bound the zero and set QINCR
	*/
	*x = small;
	/* GET-FUNCTION-VALUE */	i99999 = 1; goto L_End;
    }
    else { /* (*status > 0) */

	switch((int)i99999){
	case 1:
	    DBGprt1(" E00(i=1) ");
	    fsmall = *fx;
	    *x = big;
	    /* GET-FUNCTION-VALUE */	i99999 = 2; goto L_End;

	case 2:
	    DBGprt1(" E00(i=2) ");
	    fbig = *fx;
	    qincr = fbig > fsmall;
	    if(qincr) {
		if(fsmall > 0) {
		    *status = -1; *qleft = *qhi = 1; return;
		}
		if(fbig < 0) {
		    *status = -1; *qleft = *qhi = 0; return;
		}
	    }
	    else {

		if(fsmall < 0) {
		    *status = -1; *qleft = 1; *qhi = 0; return;
		}
		if(fbig > 0) {
		    *status = -1; *qleft = 0; *qhi = 1; return;
		}
	    }

	    *x = xsave;
	    step = fmax2(absstp,relstp*fabs(*x));
	    /* YY = F(X) - Y */
	    /* GET-FUNCTION-VALUE */	i99999 = 3; goto L_End;

	case 3: goto S90;
	case 4: goto S130;/* inmidst of `if' of S90 .. what a mess */
	case 5: goto S200;
	case 6: goto S270;
	default: return;
	}
    }

S90:
    DBGprt1(" E00(i=3) ");
    yy = *fx;
    if(yy == 0) {
	*status = 0;
	qok = 1;
	return;
    }

    qup = ( qincr && yy < 0) ||
	  (!qincr && yy > 0);
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     HANDLE CASE IN WHICH WE MUST STEP HIGHER
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
    if(!qup) goto S170;

    xlb = xsave;
    xub = fmin2(xlb+step,big);

  S120:
    /* YY = F(XUB) - Y */	*x = xub;
    /* GET-FUNCTION-VALUE */	i99999 = 4;    goto L_End;

  S130:
    DBGprt1(" E00(i=4) ");
    yy = *fx;
    qbdd = (qincr && yy >= 0) ||
	  (!qincr && yy <= 0);
    qlim = xub >= big;
    qcond = qbdd || qlim;
    if(!qcond) {
	step = stpmul*step;
	xlb = xub;
	xub = fmin2(xlb+step,big);
	goto S120;
    }

    if(qlim && !qbdd) {
	*status = -1;
	*qleft = 0;
	*qhi = !qincr;
	*x = big;
	return;
    }
    goto S240;

 S170:

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      HANDLE CASE IN WHICH WE MUST STEP LOWER
      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    xub = xsave;
    xlb = fmax2(xub-step,small);

 S190:
    /* YY = F(XLB) - Y */	*x = xlb;
    /* GET-FUNCTION-VALUE */ i99999 = 5; goto L_End;

 S200:
    DBGprt1(" E00(i=5) ");
    yy = *fx;
    qbdd = ( qincr && yy <= 0) ||
	   (!qincr && yy >= 0);
    qlim = xlb <= small;
    qcond = qbdd || qlim;
    if(!qcond) {
	step = stpmul*step;
	xub = xlb;
	xlb = fmax2(xub-step,small);
	goto S190;
    }

    if(qlim && !qbdd) {
	*status = -1;
	*qleft = 1;
	*qhi = qincr;
	*x = small;
	return;
    }

S240:
    dstzr(&xlb,&xub,&abstol,&reltol);
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
    *status = 0;

S260:
    dzror(status,x,fx,&xlo,&xhi,&qdum1,&qdum2);
    if(*status == 1){
	/* GET-FUNCTION-VALUE */ i99999 = 6; goto L_End;
    }

S270:
    DBGprt1(" E00(i=6) ");

    if(*status == 1) goto S260;

    *x = xlo;
    *status = 0;    return;

L_End:
/*
     TO GET-FUNCTION-VALUE
*/
    DBGprt2(" i99= %d ",i99999);
    *status = 1;    return;

#undef is_monotone
} /* E0000() */

void dinvr(int *status,double *x,double *fx, logical *qleft,logical *qhi)
{
/*
**********************************************************************

     void dinvr(int *status,double *x,double *fx,
	   logical *qleft,logical *qhi)

	  Double precision
	  bounds the zero of the function and invokes zror
		    Reverse Communication


			      Function


     Bounds the	   function  and  invokes  ZROR	  to perform the   zero
     finding.  STINVR  must  have   been  called  before this	routine
     in order to set its parameters.


			      Arguments


     STATUS <--> At the beginning of a zero finding problem, STATUS
		 should be set to 0 and INVR invoked.  (The value
		 of parameters other than X will be ignored on this cal

		 When INVR needs the function evaluated, it will set
		 STATUS to 1 and return.  The value of the function
		 should be set in FX and INVR again called without
		 changing any of its other parameters.

		 When INVR has finished without error, it will return
		 with STATUS 0.	 In that case X is approximately a root
		 of F(X).

		 If INVR cannot bound the function, it returns status
		 -1 and sets QLEFT and QHI.
			 INTEGER STATUS

     X <-- The value of X at which F(X) is to be evaluated.
			 DOUBLE PRECISION X

     FX --> The value of F(X) calculated when INVR returns with
	    STATUS = 1.
			 DOUBLE PRECISION FX

     QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
	  case it is .TRUE. If the stepping search terminated
	  unsucessfully at SMALL.  If it is .FALSE. the search
	  terminated unsucessfully at BIG.
		    QLEFT is LOGICAL

     QHI <-- Defined only if QMFINV returns .FALSE.  In that
	  case it is .TRUE. if F(X) > Y at the termination
	  of the search and .FALSE. if F(X) < Y at the
	  termination of the search.
		    QHI is LOGICAL

**********************************************************************
*/
    DBGprt6("N> dinvr(%d, %8g,%8g, %ld,%ld)", *status, *x, *fx, *qleft,*qhi);
    E0000(0, status,x,fx, qleft,qhi,
	  NULL,NULL, NULL,NULL, NULL,NULL,NULL);
    DBGprt6("\n -->    (%d, %8g,%8g, %ld,%ld)\n", *status, *x, *fx, *qleft,*qhi);
} /* dinvr() */

void dstinv(double *zsmall,double *zbig,double *zabsst,
	    double *zrelst,double *zstpmu,double *zabsto,
	    double *zrelto)
{
/*
**********************************************************************
      void dstinv(double *zsmall,double *zbig,double *zabsst,
	    double *zrelst,double *zstpmu,double *zabsto,
	    double *zrelto)

      Double Precision - SeT INverse finder - Reverse Communication
			      Function
     Concise Description - Given a monotone function F finds X
     such that F(X) = Y.  Uses Reverse communication -- see invr.
     This routine sets quantities needed by INVR.
	  More Precise Description of INVR -
     F must be a monotone function, the results of QMFINV are
     otherwise undefined.  QINCR must be .TRUE. if F is non-
     decreasing and .FALSE. if F is non-increasing.
     QMFINV will return .TRUE. if and only if F(SMALL) and
     F(BIG) bracket Y, i. e.,
	  QINCR is .TRUE. and F(SMALL)<=Y<=F(BIG) or
	  QINCR is .FALSE. and F(BIG)<=Y<=F(SMALL)
     if QMFINV returns .TRUE., then the X returned satisfies
     the following condition.  let
	       TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
     then if QINCR is .TRUE.,
	  F(X-TOL(X)) <= Y <= F(X+TOL(X))
     and if QINCR is .FALSE.
	  F(X-TOL(X)) >= Y >= F(X+TOL(X))
			      Arguments
     SMALL --> The left endpoint of the interval to be
	  searched for a solution.
		    SMALL is DOUBLE PRECISION
     BIG --> The right endpoint of the interval to be
	  searched for a solution.
		    BIG is DOUBLE PRECISION
     ABSSTP, RELSTP --> The initial step size in the search
	  is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
		    ABSSTP is DOUBLE PRECISION
		    RELSTP is DOUBLE PRECISION
     STPMUL --> When a step doesn't bound the zero, the step
		size is multiplied by STPMUL and another step
		taken.	A popular value is 2
		    DOUBLE PRECISION STPMUL
     ABSTOL, RELTOL --> Two numbers that determine the accuracy
	  of the solution.  See function for a precise definition.
		    ABSTOL is DOUBLE PRECISION
		    RELTOL is DOUBLE PRECISION
			      Method
     Compares F(X) with Y for the input value of X then uses QINCR
     to determine whether to step left or right to bound the
     desired x.	 the initial step size is
	  MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
     Iteratively steps right or left until it bounds X.
     At each step which doesn't bound X, the step size is doubled.
     The routine is careful never to step beyond SMALL or BIG.	If
     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
     after setting QLEFT and QHI.
     If X is successfully bounded then Algorithm R of the paper
     'Two Efficient Algorithms with Guaranteed Convergence for
     Finding a Zero of a Function' by J. C. P. Bus and
     T. J. Dekker in ACM Transactions on Mathematical
     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
     to find the zero of the function F(X)-Y. This is routine
     QRZERO.
**********************************************************************
*/
    DBGprt5("N> dstinv(%g,%g, %g,%g", *zabsst,*zabsto, *zbig,*zrelst);
    DBGprt4(" ,%g,%g, %g)\n", *zrelto,*zsmall, *zstpmu);
    E0000(1, NULL,NULL,NULL, NULL,NULL,
	  zabsst,zabsto, zbig,zrelst, zrelto,zsmall, zstpmu);
} /* dstinv() */

double dt1(double *p,double *q,double *df)
{
/*
**********************************************************************

     double dt1(double *p,double *q,double *df)
     Double precision Initalize Approximation to
	   INVerse of the cumulative T distribution


			      Function


     Returns  the  inverse   of	 the T	 distribution	function, i.e.,
     the integral from 0 to INVT of the T density is P. This is an
     initial approximation


			      Arguments


     P --> The p-value whose inverse from the T distribution is
	  desired.
		    P is DOUBLE PRECISION

     Q --> 1-P.
		    Q is DOUBLE PRECISION

     DF --> Degrees of freedom of the T distribution.
		    DF is DOUBLE PRECISION

**********************************************************************
*/
static double coef[4][5] = {
    {1,1,0,0,0},
    {3,16,5,0,0},
    {-15,17,19,3,0},
    {-945,-1920,1482,776,79}
};
static double denom[4] = {
    4,96,384,92160
};
static int ideg[4] = {
    2,3,4,5
};
static double dt1,denpow,sum,term,x,xp,xx;
static int i;
/*
     ..
     .. Executable Statements ..
*/
    x = fabs(dinvnr(p,q));
    xx = x*x;
    sum = x;
    denpow = 1;
    for(i=0; i<4; i++) {
	term = devlpl(&coef[i][0],ideg[i],xx)*x;
	denpow *= *df;
	sum += (term/(denpow*denom[i]));
    }
    if(!(*p >= 0.5)) goto S20;
    xp = sum;
    goto S30;
S20:
    xp = -sum;
S30:
    dt1 = xp;
    return dt1;
}

/* DEFINE DZROR */
static void E0001(int IENTRY,int *status,double *x,double *fx,
		  double *xlo,double *xhi,logical *qleft,
		  logical *qhi,double *zabstl,double *zreltl,
		  double *zxhi,double *zxlo)
{
#define ftol(zx) (0.5*fmax2(abstol,reltol*fabs((zx))))
  static double a,abstol,b,c,d,fa,fb,fc,fd,fda,fdb,m,mb,p,q,reltol,tol,w,xxhi,xxlo;
  static int ext,i99999;
  static logical first,qrzero;
  switch(IENTRY){
  case 0: goto DZROR;
  case 1: goto DSTZR;
  }

 DZROR:
  if(*status > 0) goto S280;
  *xlo = xxlo;
  *xhi = xxhi;
  b = *x = *xlo;
  /*
    GET-FUNCTION-VALUE
  */
  i99999 = 1;
  goto S270;
 S10:
  fb = *fx;
  *xlo = *xhi;
  a = *x = *xlo;
  /*
    GET-FUNCTION-VALUE
  */
  i99999 = 2;
  goto S270;
 S20:
  /*
    Check that F(ZXLO) < 0 < F(ZXHI)  or
    F(ZXLO) > 0 > F(ZXHI)
  */
  if(!(fb < 0)) goto S40;
  if(!(*fx < 0)) goto S30;
  *status = -1;
  *qleft = *fx < fb;
  *qhi = 0;
  return;
 S40:
 S30:
  if(!(fb > 0)) goto S60;
  if(!(*fx > 0)) goto S50;
  *status = -1;
  *qleft = *fx > fb;
  *qhi = 1;
  return;
 S60:
 S50:
  fa = *fx;
  first = 1;
 S70:
  c = a;
  fc = fa;
  ext = 0;
 S80:
  if(!(fabs(fc) < fabs(fb))) goto S100;
  if(!(c != a)) goto S90;
  d = a;
  fd = fa;
 S90:
  a = b;
  fa = fb;
  *xlo = c;
  b = *xlo;
  fb = fc;
  c = a;
  fc = fa;
 S100:
  tol = ftol(*xlo);
  m = (c+b)*.5;
  mb = m-b;
  if(!(fabs(mb) > tol)) goto S240;
  if(!(ext > 3)) goto S110;
  w = mb;
  goto S190;
 S110:
  tol = fifdsign(tol,mb);
  p = (b-a)*fb;
  if(!first) goto S120;
  q = fa-fb;
  first = 0;
  goto S130;
 S120:
  fdb = (fd-fb)/(d-b);
  fda = (fd-fa)/(d-a);
  p = fda*p;
  q = fdb*fa-fda*fb;
 S130:
  if(!(p < 0)) goto S140;
  p = -p;
  q = -q;
 S140:
  if(ext == 3) p *= 2;
  if(!(p*1 == 0 || p <= q*tol)) goto S150;
  w = tol;
  goto S180;
 S150:
  if(!(p < mb*q)) goto S160;
  w = p/q;
  goto S170;
 S160:
  w = mb;
 S190:
 S180:
 S170:
  d = a;
  fd = fa;
  a = b;
  fa = fb;
    b += w;
    *xlo = b;
    *x = *xlo;
/*
     GET-FUNCTION-VALUE
*/
    i99999 = 3;
    goto S270;
S200:
    fb = *fx;
    if(!(fc*fb >= 0)) goto S210;
    goto S70;
S210:
    if(!(w == mb)) goto S220;
    ext = 0;
    goto S230;
S220:
    ext += 1;
S230:
    goto S80;
S240:
    *xhi = c;
    qrzero = (fc >= 0 && fb <= 0) ||
	     (fc <  0 && fb >= 0);
    if(qrzero)
      *status = 0;
    else
      *status = -1;
    return;

/*----------*/
DSTZR:
/*----------*/
    xxlo = *zxlo;
    xxhi = *zxhi;
    abstol = *zabstl;
    reltol = *zreltl;
    return;
S270:
/*
     TO GET-FUNCTION-VALUE
*/
    *status = 1;
    return;
S280:
    switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S200;
      default: break;}
#undef ftol
}

void dzror(int *status,double *x,double *fx,double *xlo,
	   double *xhi,logical *qleft,logical *qhi)
/*
**********************************************************************

     void dzror(int *status,double *x,double *fx,double *xlo,
	   double *xhi,logical *qleft,logical *qhi)

     Double precision ZeRo of a function -- Reverse Communication


			      Function


     Performs the zero finding.	 STZROR must have been called before
     this routine in order to set its parameters.


			      Arguments


     STATUS <--> At the beginning of a zero finding problem, STATUS
		 should be set to 0 and ZROR invoked.  (The value
		 of other parameters will be ignored on this call.)

		 When ZROR needs the function evaluated, it will set
		 STATUS to 1 and return.  The value of the function
		 should be set in FX and ZROR again called without
		 changing any of its other parameters.

		 When ZROR has finished without error, it will return
		 with STATUS 0.	 In that case (XLO,XHI) bound the answe

		 If ZROR finds an error (which implies that F(XLO)-Y an
		 F(XHI)-Y have the same sign, it returns STATUS -1.  In
		 this case, XLO and XHI are undefined.
			 INTEGER STATUS

     X <-- The value of X at which F(X) is to be evaluated.
			 DOUBLE PRECISION X

     FX --> The value of F(X) calculated when ZROR returns with
	    STATUS = 1.
			 DOUBLE PRECISION FX

     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
	     inverval in X containing the solution below.
			 DOUBLE PRECISION XLO

     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
	     inverval in X containing the solution above.
			 DOUBLE PRECISION XHI

     QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
		at XLO.	 If it is .FALSE. the search terminated
		unsucessfully at XHI.
		    QLEFT is LOGICAL

     QHI <-- .TRUE. if F(X) > Y at the termination of the
	      search and .FALSE. if F(X) < Y at the
	      termination of the search.
		    QHI is LOGICAL

**********************************************************************
*/
{
    E0001(0,status,x,fx,xlo,xhi,qleft,qhi,NULL,NULL,NULL,NULL);
}
void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl)
/*
**********************************************************************
     void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl)
     Double precision SeT ZeRo finder - Reverse communication version
			      Function
     Sets quantities needed by ZROR.  The function of ZROR
     and the quantities set is given here.
     Concise Description - Given a function F
     find XLO such that F(XLO) = 0.
	  More Precise Description -
     Input condition. F is a double precision function of a single
     double precision argument and XLO and XHI are such that
	  F(XLO)*F(XHI)	 <=  0
     If the input condition is met, QRZERO returns .TRUE.
     and output values of XLO and XHI satisfy the following
	  F(XLO)*F(XHI)	 <= 0.
	  ABS(F(XLO)  <= ABS(F(XHI)
	  ABS(XLO-XHI)	<= TOL(X)
     where
	  TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
     If this algorithm does not find XLO and XHI satisfying
     these conditions then QRZERO returns .FALSE.  This
     implies that the input condition was not met.
			      Arguments
     XLO --> The left endpoint of the interval to be
	   searched for a solution.
		    XLO is DOUBLE PRECISION
     XHI --> The right endpoint of the interval to be
	   for a solution.
		    XHI is DOUBLE PRECISION
     ABSTOL, RELTOL --> Two numbers that determine the accuracy
		      of the solution.	See function for a
		      precise definition.
		    ABSTOL is DOUBLE PRECISION
		    RELTOL is DOUBLE PRECISION
			      Method
     Algorithm R of the paper 'Two Efficient Algorithms with
     Guaranteed Convergence for Finding a Zero of a Function'
     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
     Mathematical Software, Volume 1, no. 4 page 330
     (Dec. '75) is employed to find the zero of F(X)-Y.
**********************************************************************
*/
{
    E0001(1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,zabstl,zreltl,zxhi,zxlo);
}
double erf1(double *x)
/*
-----------------------------------------------------------------------
	     EVALUATION OF THE REAL ERROR FUNCTION
-----------------------------------------------------------------------
*/
{
static const double
    c = .564189583547756,
    a[5] = {
	.771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
	.479137145607681e-01,.128379167095513 },
    b[3] = {
	.301048631703895e-02,.538971687740286e-01,.375795757275549 },
    p[8] = {
	-1.36864857382717e-07,5.64195517478974e-01,7.21175825088309,
	4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
	4.51918953711873e+02,3.00459261020162e+02 },
    q[8] = {
	1.00000000000000,1.27827273196294e+01,7.70001529352295e+01,
	2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
	7.90950925327898e+02,3.00459260956983e+02 },
    r[5] = {
	2.10144126479064,2.62370141675169e+01,2.13688200555087e+01,
	4.65807828718470,2.82094791773523e-01 },
    s[4] = {
	9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
	1.80124575948747e+01 };

 double erf1,ax,bot,t,top,x2;
/*
     ..
     .. Executable Statements ..
*/
    ax = fabs(*x);
    if(ax > 0.5) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1;
    erf1 = *x*(top/bot);
    return erf1;
S10:
    if(ax > 4) goto S20;
    top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+
	p[7];
    bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+
	q[7];
    erf1 = 0.5+(0.5-exp(-(*x**x))*top/bot);
    if(*x < 0) return(-erf1);
    else return(erf1);
S20:
    if(ax >= 5.8) goto S30;
    x2 = *x**x;
    t = 1/x2;
    top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
    bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1;
    erf1 = (c-top/(x2*bot))/ax;
    erf1 = 0.5+(0.5-exp(-x2)*erf1);
    if(*x < 0) erf1 = -erf1;
    return erf1;
S30:
    erf1 = fifdsign(1,*x);
    return erf1;
}

double erfc1(int *ind,double *x)
/*
-----------------------------------------------------------------------
	 EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

	  ERFC1(IND,X) = ERFC(X)	    IF IND = 0
	  ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
-----------------------------------------------------------------------
*/
{
static double c = .564189583547756;
static double a[5] = {
    .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
    .479137145607681e-01,.128379167095513
};
static double b[3] = {
    .301048631703895e-02,.538971687740286e-01,.375795757275549
};
static double p[8] = {
    -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309,
    4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
    4.51918953711873e+02,3.00459261020162e+02
};
static double q[8] = {
    1.00000000000000,1.27827273196294e+01,7.70001529352295e+01,
    2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
    7.90950925327898e+02,3.00459260956983e+02
};
static double r[5] = {
    2.10144126479064,2.62370141675169e+01,2.13688200555087e+01,
    4.65807828718470,2.82094791773523e-01
};
static double s[4] = {
    9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
    1.80124575948747e+01
};
static double erfc1,ax,bot,e,t,top,w;
/*
     ..
     .. Executable Statements ..
*/
/*
		     ABS(X) <= 0.5
*/
    ax = fabs(*x);
    if(ax > 0.5) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1;
    erfc1 = 0.5+(0.5-*x*(top/bot));
    if(*ind != 0) erfc1 = exp(t)*erfc1;
    return erfc1;
S10:
/*
		  0.5 < ABS(X) <= 4
*/
    if(ax > 4) goto S20;
    top=((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[7];
    bot=((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[7];
    erfc1 = top/bot;
    goto S40;
S20:
/*
		      ABS(X) > 4
*/
    if(*x <= -5.6) goto S60;
    if(*ind != 0) goto S30;
    if(*x > 100) goto S70;
    if(*x**x > -exparg(1)) goto S70;
S30:
    t = pow(1/ *x,2);
    top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
    bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1;
    erfc1 = (c-t*top/bot)/ax;
S40:
/*
		      FINAL ASSEMBLY
*/
    if(*ind == 0) goto S50;
    if(*x < 0) erfc1 = 2*exp(*x**x)-erfc1;
    return erfc1;
S50:
    w = *x**x;
    t = w;
    e = w-t;
    erfc1 = (0.5+(0.5-e))*exp(-t)*erfc1;
    if(*x < 0) erfc1 = 2-erfc1;
    return erfc1;
S60:
/*
	     LIMIT VALUE FOR LARGE NEGATIVE X
*/
    erfc1 = 2;
    if(*ind != 0) erfc1 = 2*exp(*x**x);
    return erfc1;
S70:
/*
	     LIMIT VALUE FOR LARGE POSITIVE X
		       WHEN IND = 0
*/
    erfc1 = 0;
    return erfc1;
} /* erfc1() */

double esum(int *mu,double *x)
/*
-----------------------------------------------------------------------
		    EVALUATION OF EXP(MU + X)
-----------------------------------------------------------------------
*/
{
static double esum,w;
/*
     ..
     .. Executable Statements ..
*/
    if(*x > 0) goto S10;
    if(*mu < 0) goto S20;
    w = (double)*mu+*x;
    if(w > 0) goto S20;
    esum = exp(w);
    return esum;
S10:
    if(*mu > 0) goto S20;
    w = (double)*mu+*x;
    if(w < 0) goto S20;
    esum = exp(w);
    return esum;
S20:
    w = *mu;
    esum = exp(w)*exp(*x);
    return esum;
}
double exparg(int l)
/*
--------------------------------------------------------------------
     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
     EXP(W) CAN BE COMPUTED.

     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.

     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
--------------------------------------------------------------------
*/
{
    double lnb;
    int b; int m;

    b = ipmpar(4);

    if(b == 2)		lnb = .69314718055995;
    else if(b == 8)	lnb = 2.0794415416798;
    else if(b == 16)	lnb = 2.7725887222398;
    else		lnb = log((double)b);

    if(l == 0)	m = ipmpar(10);
    else	m = ipmpar(9)-1;

    return 0.99999*((double)m*lnb);
}

double fpser(double *a,double *b,double *x,double *eps)
/*
-----------------------------------------------------------------------

		 Evaluation of I (a,b)
				x

	  for b < min(eps,eps*a) and x <= 0.5.

-----------------------------------------------------------------------
*/
{
    double fpser,an,c,s,t,tol;
/*
	fpser := x ^ a
*/
    if(*a <= 1.e-3 * *eps) {
	fpser = 1;
    } else {
	t = *a*log(*x);
	if(t < exparg(1)) return 0;
	fpser = exp(t);
    }
/*
	Note that 1/B(a,b) = b
*/
    fpser = *b/ *a*fpser;
    tol = *eps/ *a;
    an = *a+1;
    t = *x;
    s = t/an;

    do {
	an++;
	t *= *x;
	c = t/an;
	s += c;
    } while(fabs(c) > tol);
    fpser *= (1+*a*s);
    return fpser;
} /* fpser() */

double gam1(double *a)
{
/*
     ------------------------------------------------------------------
     Computation of  1/gamma(a+1) - 1	for  -0.5 <= a <= 1.5
     ------------------------------------------------------------------
*/
    /*---------- CONSTANTS ----------*/
    static double s1 = .273076135303957;
    static double s2 = .559398236957378e-01;
    static double p[7] = {
	.577215664901533,-.409078193005776,-.230975380857675,
	.597275330452234e-01,.766968181649490e-02,-.514889771323592e-02,
	.589597428611429e-03
    };
    static double q[5] = {
	.100000000000000e+01,.427569613095214,.158451672430138,
	.261132021441447e-01,.423244297896961e-02
    };
    static double r[9] = {
	-.422784335098468,-.771330383816272,-.244757765222226,
	.118378989872749,.930357293360349e-03,-.118290993445146e-01,
	.223047661158249e-02,.266505979058923e-03,-.132674909766242e-03
    };

    double bot,d,t,top,w;

    t = *a;
    d = *a-0.5;
    if(d > 0) t = d-0.5;
    if(t == 0) return 0;
    if(t > 0) {
	top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
	bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1;
	w = top/bot;
	if(d <= 0)
	    return *a*w;
	else
	    return t / *a*(w-0.5-0.5);
    }
    else { /* t < 0 */
	top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+
	       r[1])*t+ r[0];
	bot = (s2*t+s1)*t+1;
	w = top/bot;
	if(d <= 0)
	    return *a*(w+0.5+0.5);
	else
	    return t*w / *a;
    }
} /* gam1() */

/* ===========> gaminv()  moved to separate file,
 * ===========	~~~~~~~~  since using NSWC version
 */

double gamln(double *a)
/*
-----------------------------------------------------------------------
	    EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
-----------------------------------------------------------------------
     written by Alfred H. Morris, Jr.
	Naval Surface Warfare Center (NSWC)
	Dahlgren, Virginia
*/
{
static double d = .41893853320467274;/* D = 0.5*(LN(2*PI) - 1) */
static double t,w;
static int i,n;
static double T1;

    if(*a <= 0.8)
	return gamln1(a)-log(*a);

    if(*a <= 2.25) {
	t = *a-0.5-0.5;
	return gamln1(&t);
    }
    if(*a < 10) {
	n = (*a - 1.25);
	t = *a;
	w = 1;
	for(i=0; i<n; i++) {
	    t -= 1;
	    w = t*w;
	}
	T1 = t-1;
	return gamln1(&T1)+log(w);
    }
    /* a >= 10 -- use DEL() as	algdiv() & bcorr() do : */
    t = pow(1/ *a,2);
    w = (((((LG_c5*t+LG_c4)*t+LG_c3)*t+LG_c2)*t+LG_c1)*t+LG_c0)/ *a;
    return (d+w+(*a-0.5)*(log(*a)-1));
} /* gamln() */

double gamln1(double *a)
/*
-----------------------------------------------------------------------
     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25
-----------------------------------------------------------------------
*/
{
    static double p0 = .577215664901533;
    static double p1 = .844203922187225;
    static double p2 = -.168860593646662;
    static double p3 = -.780427615533591;
    static double p4 = -.402055799310489;
    static double p5 = -.673562214325671e-01;
    static double p6 = -.271935708322958e-02;
    static double q1 = .288743195473681e+01;
    static double q2 = .312755088914843e+01;
    static double q3 = .156875193295039e+01;
    static double q4 = .361951990101499;
    static double q5 = .325038868253937e-01;
    static double q6 = .667465618796164e-03;
    static double r0 = .422784335098467;
    static double r1 = .848044614534529;
    static double r2 = .565221050691933;
    static double r3 = .156513060486551;
    static double r4 = .170502484022650e-01;
    static double r5 = .497958207639485e-03;
    static double s1 = .124313399877507e+01;
    static double s2 = .548042109832463;
    static double s3 = .101552187439830;
    static double s4 = .713309612391000e-02;
    static double s5 = .116165475989616e-03;

    static double w,x;
/*
     ..
     .. Executable Statements ..
*/
    if(*a < 0.6) {
	w = ((((((p6**a+p5)**a+p4)**a+p3)**a+p2)**a+p1)**a+p0)/
	    ((((((q6**a+q5)**a+q4)**a+q3)**a+q2)**a+q1)**a+1);
	return -(*a*w);
    }
    else{
	x = *a-0.5-0.5;
	w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/
	    (((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x +1);
	return x*w;
    }
} /* gamln1 */

double Xgamm(double *a) /* not gamma() [namespace!] */
{
/*
-----------------------------------------------------------------------

	 Evaluation of the Gamma function for real arguments

			   -----------

     Gamma(A) := 0   when the gamma function cannot be computed.

-----------------------------------------------------------------------
     written by Alfred H. Morris, Jr.
	Naval Surface Warfare Center (NSWC)
	Dahlgren, Virginia
-----------------------------------------------------------------------
*/
    static double d = .41893853320467274178;
    static double pi = 3.1415926535898;
    static double r1 = .820756370353826e-03;
    static double r2 = -.595156336428591e-03;
    static double r3 = .793650663183693e-03;
    static double r4 = -.277777777770481e-02;
    static double r5 = .833333333333333e-01;
    static double p[7] = {
	.539637273585445e-03,.261939260042690e-02,.204493667594920e-01,
	.730981088720487e-01,.279648642639792,.553413866010467,1
    };
    static double q[7] = {
	-.832979206704073e-03,.470059485860584e-02,.225211131035340e-01,
	-.170458969313360,-.567902761974940e-01,.113062953091122e+01,1
    };


    double Xgamm,bot,g,lnx,s=0,t,top,w,x,z;
    int i,j,m,n;

    x = *a;
    if(fabs(*a) < 15) {
	/* ---------------------------------------------------------------
	   evaluation of gamma(a) for abs(a) < 15
	   --------------------------------------------------------------- */
	t = 1;
	m = fifidint(*a)-1;
	/*
	  let t be the product of a-j when a >= 2
	*/
	if(m >= 0) {
	    for(j=0; j<m; j++) {
		x -= 1;
		t *= x;
	    }
	    x -= 1;
	}
	else {
	    /* t := the product of a+j when a < 1 */
	    t = *a;
	    if(*a <= 0) {
		m = -m-1;
		for(j=0; j<m; j++) {
		    x += 1;
		    t *= x;
		}
		x += (0.5+0.5);
		t = x*t;
		if(t == 0) return 0.;
	    }
	    /*
	      the following code checks if 1/t can overflow.
	      this code may be omitted if desired.
	    */
	    if(fabs(t) < 1e-30) {
		if(fabs(t)*spmpar(3) <= 1.0001) return 0.;
		return 1/t;
	    }
	}

	/*
	  GAMMA(1 + X) FOR  0 <= X < 1 --- rational approx (p,q) ---
	*/
	top = p[0];
	bot = q[0];
	for(i=1; i<7; i++) {
	    top = p[i]+ x*top;
	    bot = q[i]+ x*bot;
	}
	Xgamm = top/bot;

	if(*a >= 1)
	    Xgamm *= t;
	else
	    Xgamm /= t;
	return Xgamm;
    }
    else {
	/*
	  -----------------------------------------------------------
	  EVALUATION OF GAMMA(A) FOR ABS(A) >= 15
	  -----------------------------------------------------------
	*/
	if(fabs(*a) >= 1e3) return 0.;
	if(*a <= 0) {
	    x = -*a;
	    n = (long)(x);
	    t = x-(double)n;
	    if(t > 0.9) t = 1-t;
	    s = sin(pi*t)/pi;
	    if(n%2 == 0) s = -s;
	    if(s == 0) return 0.;
	}

	/*
	  COMPUTE THE MODIFIED ASYMPTOTIC SUM
	*/
	t = 1/(x*x);
	g = ((((r1*t+r2)*t+r3)*t+r4)*t+r5)/x;
	/*
	  ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
	  BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
	*/
	lnx = log(x);
	/*
	  FINAL ASSEMBLY
	*/
	z = x;
	g = d+g+(z-0.5)*(lnx-1e0);
	w = g;
	t = g-w;
	if(w > 0.99999*exparg(0)) return 0.;
	Xgamm = exp(w)*(1+t);
	if(*a < 0) Xgamm = 1/(Xgamm*s)/x;
	return Xgamm;
    }
} /* Xgamm */

void grat1(double *a,double *x,double *r,double *p,double *q, double *eps)
{
/*
-----------------------------------------------------------------------
	EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
		      P(A,X) AND Q(A,X)
     IT IS ASSUMED THAT A <= 1.	 EPS IS THE TOLERANCE TO BE USED.
     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X^A/GAMMA(A).
-----------------------------------------------------------------------
*/

static int K0 = 0;
    double a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,t,tol,w,z,T1;

    if(*a**x == 0) goto S120;
    if(*a == 0.5) {
	T1 = sqrt(*x);
	if(*x < 0.25) {
	    *p = erf1(&T1);
	    *q = 0.5+(0.5-*p);
	} else {
	    *q = erfc1(&K0,&T1);
	    *p = 0.5+(0.5-*q);
	}
	return;
    }
    else if(*x < 1.1) {
/*
  TAYLOR SERIES FOR P(A,X)/X^A
*/
	an = 3;
	c = *x;
	sum = *x/(*a+3);
	tol = 0.1**eps/(*a+1);
	do {
	    an += 1;
	    c = -(c*(*x/an));
	    t = c/(*a+an);
	    sum += t;
	} while(fabs(t) > tol);

	j = *a**x*((sum/6-0.5/(*a+2))**x+1/(*a+1));
	z = *a*log(*x);
	h = gam1(a);
	g = 1+h;
	if(*x >= 0.25) {
	    if(*a < *x/2.59) goto S50;
	}
	else {
	    if(z > -.13394) goto S50;
	}

	w = exp(z);
	*p = w*g*(0.5+(0.5-j));
	*q = 0.5+(0.5-*p);
	return;
      S50:
	l = cdf_rexp(&z);
	w = 0.5+(0.5+l);
	*q = (w*j-l)*g-h;
	if(*q < 0) goto L_ret_1_0;
	*p = 0.5+(0.5-*q);
	return;
    }

    else { /* x >= 1.1 */
/*
  CONTINUED FRACTION EXPANSION
*/
	a2nm1 = a2n = 1;
	b2nm1 = *x;
	b2n = *x+(1-*a);
	c = 1;

	do {
	    a2nm1 = *x*a2n + c*a2nm1;
	    b2nm1 = *x*b2n + c*b2nm1;
	    am0 = a2nm1/b2nm1;
	    c += 1;
	    cma = c-*a;
	    a2n = a2nm1 + cma*a2n;
	    b2n = b2nm1 + cma*b2n;
	    an0 = a2n/b2n;
	} while(fabs(an0-am0) >= *eps*an0);
	*q = *r*an0;
	*p = 0.5+(0.5-*q);
	return;
    }
/*
		SPECIAL CASES
*/
S120:
    if(*x <= *a) goto L_ret_0_1;
    goto L_ret_1_0;

L_ret_0_1:	  *p = 0;    *q = 1;	  return;
L_ret_1_0:	  *p = 1;    *q = 0;	  return;

} /* grat1() */

/* The NSWC library (1993) has an improved version of gratio()
 *     ------------
 * --> /usr/local/lib.src/NSWC/nswc/
 * and /usr/local/lib.src/NSWC/C-src/ gratio-1.c
 *
 * I have converted them via f2c & f2c-clean and by hand:
 *			Martin Maechler, Statistik, ETH Zurich, May 1999
 *
 * The `older' version in dcdflib (1.1; Dec.1997) is now in
 * ./dcdf-gratio.c
 */
void gratio(double *a, double *x,
	    double *ans, double *qans, int *ind)
{
/*
 ----------------------------------------------------------------------
	Evaluation of the INCOMPLETE GAMMA RATIO functions

		      P(A,X) AND Q(A,X)

			----------

	It is assumed that a and x are nonnegative, where a and x
	are not both 0.

	ans and qans are variables.
	gratio assigns ans the value p(a,x) and qans the value q(a,x).

	ind may be any integer.
	if ind = 0 then the user is requesting as much accuracy as
		possible (up to 14 significant digits).	 Otherwise,
	if ind = 1 then accuracy is requested to within 1 unit of the
		6-th significant digit, and
	if ind != 0,1 then accuracy
		is requested to within 1 unit of the 3rd significant digit.

	error return ...
	   ans is assigned the value 2 when a or x is negative,
	when a*x = 0, or when p(a,x) and q(a,x) are indeterminant.
	p(a,x) and q(a,x) are computationally indeterminant when
	x is exceedingly close to a and a is extremely large.
 ----------------------------------------------------------------------
     written by Alfred H. Morris, Jr.
	Naval Surface Warfare Center (NSWC)
	Dahlgren, Virginia
 -------------------------
     revised ... dec 1991

	     -------------------------------
 >>>>>>>>>>> PART of PUBLIC LIBRARY nswclib <<<<<<<<<<<
	     -------------------------------
*/

/* --------- CONSTANTS -------------------- */

static const double
    d10 = -.00185185185185185,
    d20 =  .00413359788359788,
    d30 =  6.49434156378601e-4,
    d40 = -8.61888290916712e-4,
    d50 = -3.36798553366358e-4,
    d60 =  5.31307936463992e-4,
    d70 =  3.44367606892378e-4,
    d80 = -6.52623918595309e-4,

    /* alog10= 2.30258509299405		=  M_LN10 = ln(10) */
    /* rt2pin= .398942280401432677939946 = M_1_SQRT_2PI = 1/sqrt(2 pi) */
    /* rtpi  = 1.772453850905516 	=  M_SQRT_PI = sqrt(pi) */
    third = 1/3.,

    /* These are chosen depending on ind = 0, 1, 2 : */
    acc0[3] = { 5e-15, 5e-7, 5e-4 },
    e0  [3] = { 2.5e-4, .025, .14 },
    big [3] = { 25., 14., 10. },
    x0  [3] = { 31., 17.,  9.7 },

    a1[4] = {
	-3.9878392437077e-6,-5.87926036018402e-4,
	-.0049168713172692,-.00185185185184291 },
    b1[4] = {
	.00386325038602125,.0506042559238939,
	.283344278023803,.780110511677243 },
    a2[2] = { 6.69564126155663e-4,.00413359788442192 },
    b2[5] = {
	-4.21924263980656e-4,.00650837693041777,
	.0682034997401259,.339173452092224,.810647620703045 },
    a3[2] = { 8.10586158563431e-4,6.4943415761977e-4 },
    b3[5] = {
	-6.3227658735212e-4,.00905375887385478,
	.0906610359762969,.406288930253881,.894800593794972 },
    a4[2] = { -1.05014537920131e-4,-8.61888301199388e-4 },
    b4[4] = {
	.0322609381345173,.17829577356297,
	.591353097931237,1.03151890792185 },
    a5[2] = { -4.35211415445014e-4,-3.36806989710598e-4 },
    b5[3] = { .178716720452422,.600380376956324, 1.08515217314415 },
    a6[2] = { -1.82503596367782e-4,5.31279816209452e-4 },
    b6[2] = { .345608222411837,.770341682526774 },
    a7[2] = { 4.43219646726422e-4,3.44430064306926e-4 },
    b7[2] = { .821824741357866,1.15029088777769 },
    a8[2] = { 8.78371203603888e-4,-6.86013280418038e-4 },
    d0[6] = {
	.0833333333333333,-.0148148148148148,
	.00115740740740741,3.52733686067019e-4,-1.78755144032922e-4,
	3.91926317852244e-5 },
    d1[4] = {
	-.00347222222222222,.00264550264550265,
	-9.9022633744856e-4,2.05761316872428e-4 },
    d2[2] = { -.00268132716049383,7.71604938271605e-4 },
    d3[2] = { 2.29472093621399e-4,-4.69189494395256e-4 },
    d4[1] = { 7.84039221720067e-4 },
    d5[1] = { -6.97281375836586e-5 },
    d6[1] = { -5.92166437353694e-4 },

    a0[4] = {
	-.00231272501940775, -.033537852002422,
	-.15984014344399, -.333333333333333 },
    b0[6] = {
	6.33763414209504e-7, -9.39001940478355e-6, .00239521354917408,
	.0376245718289389, .238549219145773, .729520430331981 };


/* --------- Variables -------------------- */

    double wk[20];/* Work Array */

    double a2nm1, b2nm1, twoa, c, e, g,h, j,l, r,s, t = 0. /*Wall*/;
    double u, w, y, z, c0, c1, c2, c3, c4, c5, c6, c7, c8;
    double a2n, b2n, acc, amn, apn, rta, tol, sum, rtx, T6;

    long i, m, n, iop;

/* ------------------------- */

    if(*a < 0 || *x < 0 ||(*a == 0 && *x == 0)) goto L_err_ret;
    if(*a**x == 0) goto L_0__1;

    iop = *ind;
    if(iop != 0 && iop != 1) iop = 2;/* CASES 0, 1 & 2 */
    /*
      E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
      FLOATING POINT NUMBER FOR WHICH 1 + E > 1 .
    */
    e = spmpar(1);
    acc = fmax2(acc0[iop],e);

    DBGprt3("gratio(a=%4g,x=%4g): ", *a, *x);


    if (*a < 1) {

	if (*a == 0.5) { /* L320: */
	    T6 = sqrt(*x);
	    if (*x >= 0.25) { /* L321: */
		*qans = erfc1(&I0, &T6);
		*ans = (.5 - *qans) + .5;
	    }
	    else {
		*ans = erf1(&T6);
		*qans = (.5 - *ans) + .5;
	    }
	    return;
	}
	/* 0 < a < 1 (a != 1/2) ;  x > 0 */

	if (*x < 1.1) {/* L110:		TAYLOR SERIES FOR P(A,X)/X^A */

	    l = 3;
	    c = *x;
	    sum = *x / (*a + 3);
	    tol = 3 * acc / (*a + 1);

	    do {
		l++;
		c = -c * (*x / l);
		t = c / (*a + l);
		sum += t;
	    } while (fabs(t) > tol);
	    DBGprt3(" x < 1.1 -- Taylor after do: l= %3g, sum= %g\n", l, sum);

	    j = *a * *x * ((sum / 6 - 0.5 / (*a + 2)) * *x + 1/(*a + 1));

	    z = *a * log(*x);
	    h = gam1(a);
	    g = h + 1;

	    if((*x < 0.25 && z > -.13394) || *a < *x/2.59) { /* L135: */
		l = cdf_rexp(&z);
		w = 5 + (l + .5);
		*qans = (w * j - l) * g - h;
		if (*qans < 0) goto L_ret_1_0;
		*ans = (.5 - *qans) + .5;
	    }
	    else {
		w = exp(z);
		*ans = w * g * (.5 - j + .5);
		*qans = .5 - *ans + .5;
	    }
	    return;
	}
	/* else  x >= 1.1 ; a < 1 */
	r = rcomp(a, x);/* = exp(-x)* x^a / Gamma(a) */
	DBGprt2(" x >= 1.1 (a < 1): rcomp()=r= %8g,", r);
	if (r == 0.) goto L_ret_1_0;

	goto L_cont_frac;

    }
    else {/* a >= 1 */

	DBGprt1(" a >= 1:");
	if (*a < big[iop]) {
	    DBGprt2(" 1 <= a < big=%5g ", big[iop]);
	    if (*a <= *x  &&  *x < x0[iop]) {
		twoa = *a + *a;
		m = (long) twoa;
		if (twoa == (double)m) {
		    /*
		      FINITE SUMS FOR Q WHEN  A >= 1  AND  2*A  IS AN INTEGER
		    */
		    DBGprt2(" 2*A = integer; a=%4g\n", *a);
		    i = m / 2;
		    if(*a == (double)i) { /* L140 */
			sum = exp(-*x);
			t = sum;
			n = 1;
			c = 0;
		    }
		    else { /* L150: */
			rtx = sqrt(*x);
			sum = erfc1(&I0, &rtx);
			t = exp(-*x)/(M_SQRT_PI*rtx);
			n = 0;
			c = -0.5;
		    }
		    while(n != i) {
			n++;
			c++;
			t *= *x/c;
			sum += t;
		    }
		    *qans = sum;
		    *ans = .5 + (0.5-*qans);
		    return;
		}
	    }
	}
	else { /* *a >= big[iop] > 1 (L20) */

	    l = *x / *a;
	    DBGprt3(" a >= big=%5g: l= %8g ", big[iop], l);
	    if (l == 0) goto L_ret_0_1;
	    s = 0.5+(0.5-l);
	    z = rlog(&l);
	    if (z >= 700 / *a) { /* L330: */
		if (fabs(s) <= 2*e)
		    goto L_err_ret;
		else
		    goto L_0__1;
	    }
	    y = *a * z;
	    rta = sqrt(*a);
	    if (fabs(s) <= e0[iop] / rta) goto L_L1Temme;
	    if (fabs(s) <= .4)		goto L_TEMME;

	}
	/* L30:
	   twoa != m   ||  *a > *x  || *x >= x0[iop] */
	r = rcomp(a, x);
	DBGprt2(" L30: rcomp()=r= %8g,", r);
	if (r == 0.)	    goto L_0__1;
	if (*x <= fmax2(*a,M_LN10)) { /* L50: */
	    /*
	      TAYLOR SERIES FOR P/R
	    */
	    apn = *a + 1.;
	    t = *x / apn;
	    wk[0] = t;
	    for (n = 1; n < 20; n++) {
		apn += 1.;
		t *= (*x / apn);
		if (t <= 1e-3) break;
		wk[n] = t;
	    }
	    sum = t;
	    tol = .5 * acc;
	    do {
		apn += 1.;
		t *= *x / apn;
		sum += t;
	    } while (t > tol);
	    DBGprt3("L50: t= %9g,sum= %9g,", t,sum);
	    for(; n > 0; ) {
		sum += wk[--n];
		/*DBGprt2("%d,",n);*/
	    }
	    DBGprt2(" sum= %9g\n", sum);
	    *ans = r / *a * (sum + 1.);
	    *qans = .5+(.5-*ans);
	    return;
	}
	if (*x >= x0[iop]) {/* x > max(a, x0, ln(10));  a >= 1  (L80)

			       ASYMPTOTIC EXPANSION */
	    amn = *a - 1;
	    t = amn / *x;
	    wk[0] = t;
	    for(n=1; n < 20; n++) {
		amn--;
		t *= (amn/ *x);
		if(fabs(t) <= 1e-3) break;
		wk[n] = t;
	    }
	    sum = t;
	    DBGprt3(" Asympt t0=%g, t=%g\n", amn/ *x, t);
	    while(fabs(t) > acc) {
		amn--;
		t *= (amn/ *x);
		sum += t;
	    }

	    for(; n > 0; ) {
		sum += wk[--n];
	    }
	    *qans = r/ *x*(1+sum);
	    *ans = .5+ (.5-*qans);
	    return;
	}
    }

L_cont_frac:/*		Continued Fraction Expansion */

    DBGprt1("gratio: Cont.Frac.");
    tol = fmax2(8 * e, 4 * acc);/* revised (NSWC differs from dcdflib) */
    a2nm1 = a2n = 1;
    b2nm1 = *x;
    b2n = *x + (1. - *a);
    c = 1;
    do {
	a2nm1 = *x * a2n + c * a2nm1;
	b2nm1 = *x * b2n + c * b2nm1;
	c++;
	t = c - *a;
	a2n = a2nm1 + t * a2n;
	b2n = b2nm1 + t * b2n;

	/* Re-scale */
	a2nm1 /= b2n;
	b2nm1 /= b2n;
	a2n   /= b2n;
	b2n = 1.;
    } while(fabs(a2n - a2nm1 / b2nm1) >= tol * a2n);

    DBGprt3("-> %4g terms a2n=%11g\n", c, a2n);
    *qans = r * a2n;/* where r = rcomp(a, x) = exp(-x)* x^a / Gamma(a) */
    *ans = .5 - *qans + .5;
    return;

L_TEMME:
    DBGprt2("Gen.Temme Expansion: s= %g\n",s);
    if (fabs(s) <= 2*e  &&  *a * e * e > .00328) goto L_err_ret;

    c = exp(-y);
    T6 = sqrt(y);
    w = erfc1(&I1, &T6) * .5;
    u = 1. / *a;
    z = sqrt(z + z);
    if (l < 1) z = -z;

    switch(iop) {
    case 0: /* L210: */
	if (fabs(s) <= .001) goto L_Temme_case_0;

/*            USING THE MINIMAX APPROXIMATIONS */

	c0 = (((a0[0] * z + a0[1]) * z + a0[2]) * z + a0[3]) /
	    ((((((b0[0] * z + b0[1]) * z + b0[2]) * z + b0[3]) * z + b0[4]) * z
	      + b0[5]) * z + 1.);
	c1 = (((a1[0] * z + a1[1]) * z + a1[2]) * z + a1[3]) /
	    ((((b1[0] * z + b1[1]) * z + b1[2]) * z + b1[3]) * z + 1.);
	c2 = (a2[0] * z + a2[1]) /
	    (((((b2[0]* z + b2[1])* z + b2[2])* z + b2[3])* z + b2[4])* z + 1.);
	c3 = (a3[0] * z + a3[1]) /
	    (((((b3[0]* z + b3[1])* z + b3[2])* z + b3[3])* z + b3[4])* z + 1.);
	c4 = (a4[0] * z + a4[1]) /
	    ((((b4[0] * z + b4[1]) * z + b4[2]) * z + b4[3]) * z + 1.);
	c5 = (a5[0] * z + a5[1]) / (((b5[0] * z + b5[1]) * z + b5[2]) * z + 1.);
	c6 = (a6[0] * z + a6[1]) / ((b6[0] * z + b6[1]) * z + 1.);
	c7 = (a7[0] * z + a7[1]) / ((b7[0] * z + b7[1]) * z + 1.);
	c8 = a8[0] * z + a8[1];
	t = (((((((c8 * u + c7) * u + c6) * u + c5) * u + c4) * u + c3)
	      * u + c2) * u + c1) * u + c0;
	break;

    case 1: /* L220: */
/*                    TEMME EXPANSION */

	c0 = (((((d0[5] * z + d0[4]) * z + d0[3]) * z + d0[2]) * z + d0[1]
	    ) * z + d0[0]) * z - third;
	c1 = (((d1[3] * z + d1[2]) * z + d1[1]) * z + d1[0]) * z + d10;
	c2 = d2[0] * z + d20;
	t = (c2 * u + c1) * u + c0;

	break;

    case 2: /* L230: */
	t = ((d0[2] * z + d0[1]) * z + d0[0]) * z - third;
    }

L_End_Temme:
    if (l < 1.) { /* L241 */
	*ans = c * (w - M_1_SQRT_2PI * t / rta);
	*qans = .5 - *ans + .5;
    }
    else {
	*qans = c * (w + M_1_SQRT_2PI * t / rta);
	*ans = .5 - *qans + .5;
    }
    return;

/*               TEMME EXPANSION FOR L = 1 */

L_L1Temme:
    DBGprt4("Temme Expansion (L~=1): L= %g, y= %g, z= %g;", l, y, z);
    if(*a*e*e > 3.28e-3) goto L_err_ret;

    c = 0.5+(0.5-y);
    w = (0.5 - sqrt(y)*(0.5+(0.5 - y/3))/M_SQRT_PI)/c;
    u = 1. / *a;
    z = sqrt(z + z);
    if(l < 1) z = -z;

    switch(iop) {
    case 0:
    L_Temme_case_0:
	c0 = ((d0[2] * z + d0[1]) * z + d0[0]) * z - third;
	c1 = ((d1[2] * z + d1[1]) * z + d1[0]) * z + d10;
	c2 = (d2[1] * z + d2[0]) * z + d20;
	c3 = (d3[1] * z + d3[0]) * z + d30;
	c4 = d4[0] * z + d40;
	c5 = d5[0] * z + d50;
	c6 = d6[0] * z + d60;
	t = (((((((d80 * u + d70) * u + c6) * u + c5) * u + c4) * u + c3)
	      * u + c2) * u + c1) * u + c0;
	break;

    case 1: /* L270: */
	c0 = (d0[1] * z + d0[0]) * z - third;
	c1 = d1[0] * z + d10;
	t = (d20 * u + c1) * u + c0;
	break;

    case 2: /* L280: */
	t = d0[0] * z - third;
    }
    goto L_End_Temme;

/*
		     Some Special Cases
*/
 L_0__1:

    if(*x <= *a) {
 L_ret_0_1: DBGprt1("ret_0_1\n"); *ans = 0;   *qans = 1;	 return;
    }
    else {
 L_ret_1_0: DBGprt1("ret_1_0\n"); *ans = 1;   *qans = 0;	 return;
    }

 L_err_ret: DBGprt1("ERROR_ret\n"); 	*ans = 2; return;

} /* gratio */


double gsumln(double *a,double *b)
/*
-----------------------------------------------------------------------
	  EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
	  FOR 1 <= A <= 2  AND	1 <= B <= 2
-----------------------------------------------------------------------
*/
{
    double x,T1;

    x = *a+*b-2e0;
    if(x <= 0.25) {
	T1 = 1+x;
	return gamln1(&T1);
    }
    else if(x <= 1.25) {
	return gamln1(&x)+alnrel(&x);
    }
    else {
	T1 = x-1;
	return gamln1(&T1)+log(x*(1+x));
    }
} /* gsumln() */


double psi(double *xx)
{
/*
---------------------------------------------------------------------

		 Evaluation of the digamma function

			   -----------

     psi(xx) is assigned the value 0 when the digamma function cannot
     be computed.

     The main computation involves evaluation of rational chebyshev
     approximations published in
     Math. Comp. 27(1973), 123-127  by	Cody, Strecok and Thacher.

---------------------------------------------------------------------
     psi was written at Argonne National Laboratory for the FUNPACK
     package of special function subroutines.
     psi was modified by     A.H. Morris (NSWC).
---------------------------------------------------------------------
*/
static double dx0 = 1.461632144968362341262659542325721325;
static double piov4 = .785398163397448;
static double p1[7] = {
    .895385022981970e-02,.477762828042627e+01,.142441585084029e+03,
    .118645200713425e+04,.363351846806499e+04,.413810161269013e+04,
    .130560269827897e+04
};
static double p2[4] = {
    -.212940445131011e+01,-.701677227766759e+01,-.448616543918019e+01,
    -.648157123766197
};
static double q1[6] = {
    .448452573429826e+02,.520752771467162e+03,.221000799247830e+04,
    .364127349079381e+04,.190831076596300e+04,.691091682714533e-05
};
static double q2[4] = {
    .322703493791143e+02,.892920700481861e+02,.546117738103215e+02,
    .777788548522962e+01
};

    int i,m,n,nq;
    double aug,den,sgn,upper,w,x,xmx0,z;
    double xmax1,xsmall;

/*
---------------------------------------------------------------------
     machine dependent constants ...
	xmax1  = the smallest positive floating point constant
		 with entirely integer representation.	Also used
		 as negative of lower bound on acceptable negative
		 arguments and as the positive argument beyond which
		 psi may be represented as alog(x).
	xsmall = absolute argument below which pi*cotan(pi*x)
		 may be represented by 1/x.
---------------------------------------------------------------------
*/
    xmax1 = ipmpar(3);
    xmax1 = fmin2(xmax1,1/spmpar(1));
    xsmall = 1e-9;

    x = *xx;
    aug = 0;
    if(x < 0.5) {
/*
  ---------------------------------------------------------------------
  X < 0.5,  USE REFLECTION FORMULA
  PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
  ---------------------------------------------------------------------
*/
	if(fabs(x) <= xsmall) {
	    if(x == 0) goto Err;
	    /*
	      ---------------------------------------------------------------
	      0 < |x| <= xsmall.  Use 1/x as substitute for  PI*COTAN(PI*X)
	      ---------------------------------------------------------------
	    */
	    aug = -(1/x);
	}
	else { /* |x| > xsmall */

	    /*
	      ----------------------------------
	      REDUCTION OF ARGUMENT FOR COTAN
	      ----------------------------------
	    */
	    w = -x;
	    sgn = piov4;
	    if(w <= 0) {
		w = -w;
		sgn = -sgn;
	    }
	    /*
	      -------------------------------------
	      make an error exit if x <= -xmax1
	      -------------------------------------
	    */
	    if(w >= xmax1) goto Err;
	    nq = fifidint(w);
	    w -= (double)nq;
	    nq = fifidint(w*4);
	    w = 4*(w-(double)nq*.25);
	    /*
	      -------------------------------------------------------------
	      w is now related to the fractional part of	 4 * x.
	      adjust argument to correspond to values in first
	      quadrant and determine sign
	      -------------------------------------------------------------
	    */
	    n = nq/2;
	    if(n+n != nq) w = 1-w;
	    z = piov4*w;
	    m = n/2;
	    if(m+m != n) sgn = -sgn;
	    /*
	      ------------------------------------------------
	      Determine final value for		-pi*cotan(pi*x)
	      ------------------------------------------------
	    */
	    n = (nq+1)/2;
	    m = n/2;
	    m += m;
	    if(m == n) {
		/*
		  check for singularity
		*/
		if(z == 0) goto Err;
		/*
		  -----------------------------------------------
		  use cos/sin as a substitute for cotan, and
		  sin/cos as a substitute for tan
		  -----------------------------------------------
		*/
		aug = sgn*(cos(z)/sin(z)*4);
	    }
	    else { /* m != n */

		aug = sgn*(sin(z)/cos(z)*4);
	    }
	}

	x = 1-x;

    } /* end if (x < .5) */

    if(x <= 3) {
	/*
	  --------------------------------------
	  0.5 <= X <= 3
	  --------------------------------------
	*/
	den = x;
	upper = p1[0]*x;
	for(i=1; i<=5; i++) {
	    den = (den+q1[i-1])*x;
	    upper = (upper+p1[i+1-1])*x;
	}
	den = (upper+p1[6])/(den+q1[5]);
	xmx0 = x-dx0;
	return den*xmx0+aug;
    }

    if(x < xmax1) {
	/*
	  ----------------
	  3 < X < XMAX1
	  ----------------
	*/
	w = 1/(x*x);
	den = w;
	upper = p2[0]*w;
	for(i=1; i<=3; i++) {
	    den = (den+q2[i-1])*w;
	    upper = (upper+p2[i+1-1])*w;
	}
	aug = upper/(den+q2[3])-0.5/x+aug;
    }
    /*
      -----------------------------
      IF X >= XMAX1, PSI = LN(X)
      -----------------------------
    */

    return aug+log(x);

/*
---------------------------------------------------------------------
     ERROR RETURN
---------------------------------------------------------------------
*/
Err:
    return 0;
} /* psi() */


double rcomp(double *a, double *x)
{
/* -----------------------------------------------------------------------
		Evaluation of  EXP(-X)* X^A / GAMMA(A)
 -----------------------------------------------------------------------
*/
    double t, t1, u;

    if (*x == 0.) return(0.);

    if (*a < 20.) {
	t = *a * log(*x) - *x;
	if (t < exparg(1))
	    return(0.);/* underflow */

	if (*a < 1.)
	    return(*a * exp(t) * (gam1(a) + 1.));
	else
	    return(exp(t) / Xgamm(a));
    }
    else { /* a >= 20 */
	u = *x / *a;
	if (u == 0.) {
	    return(0.);
	}
	t = 1. / (*a * *a);
	t1 = (((t * .75 - 1.) * t + 3.5) * t - 105.) / (*a * 1260.);
	t1 -= *a * rlog(&u);
	if (t1 < exparg(1))
	    return(0.);/* underflow */
	else
	    return(M_1_SQRT_2PI * sqrt(*a) * exp(t1));
    }
} /* rcomp() */

double cdf_rexp(double *x)
/*
-----------------------------------------------------------------------
	    Evaluation of  EXP(x) - 1
-----------------------------------------------------------------------
*/
{
static double p1 = .914041914819518e-09;
static double p2 = .238082361044469e-01;
static double q1 = -.499999999085958;
static double q2 = .107141568980644;
static double q3 = -.119041179760821e-01;
static double q4 = .595130811860248e-03;

    double w;

    if(fabs(*x) <= 0.15) {
	return *x*(((p2**x+p1)**x+1)/((((q4**x+q3)**x+q2)**x+q1)**x+1));
    }
    else {
	w = exp(*x);
	if(*x <= 0)
	    return w-0.5-0.5;
	else
	    return w*(0.5+(0.5-1/w));
    }
} /* cdf_rexp() */

double rlog(double *x)
/*
	----------------------------
	Evaluation of  x - 1 - LN(x)
	----------------------------
*/
{
static double a = .566749439387324e-01;
static double b = .456512608815524e-01;
static double p0 = .333333333333333;
static double p1 = -.224696413112536;
static double p2 = .620886815375787e-02;
static double q1 = -.127408923933623e+01;
static double q2 = .354508718369557;

    double r,t,u,w,w1;

    if(*x < 0.61 || *x > 1.57) {
	r = *x-0.5-0.5;
	return r-log(*x);
    }
    else { /*========	0.61 <= x <= 1.57 :

	Argument Reduction */

	if(*x < 0.82) { /* S10 */
	    u = *x-0.7;
	    u /= 0.7;
	    w1 = a-u*0.3;
	}
	else if(*x > 1.18) { /* S20 */
	    u = 0.75**x-1e0;
	    w1 = b+u/3;
	}
	else {
	    u = *x-0.5-0.5;
	    w1 = 0;
	}
	/*	Series Expansion */
	r = u/(u+2);
	t = r*r;
	w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1);
	return 2*t*(1/(1-r)-r*w)+w1;
    }
} /* rlog() */

double rlog1(double *x)
/*
-----------------------------------------------------------------------
	     Evaluation of  X - LN(1 + X)
-----------------------------------------------------------------------
*/
{
static double a = .566749439387324e-01;
static double b = .456512608815524e-01;
static double p0 = .333333333333333;
static double p1 = -.224696413112536;
static double p2 = .620886815375787e-02;
static double q1 = -.127408923933623e+01;
static double q2 = .354508718369557;

    double h,r,t,w,w1;

    if(*x < -0.39 || *x > 0.57) {
	w = *x+0.5+0.5;
	return *x-log(w);
    }

/*========= -0.39 <= x <= 0.57 =======

	Argument Reduction
*/
    if(*x < -0.18) {
	h = *x+0.3;
	h /= 0.7;
	w1 = a-h*0.3;
    } else if(*x > 0.18) {
	h = 0.75**x-0.25;
	w1 = b+h/3;
    } else {
	h = *x;
	w1 = 0;
    }

/*
	Series Expansion
*/
    r = h/(h+2);
    t = r*r;
    w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1);
    return 2*t*(1/(1-r)-r*w)+w1;
}

double spmpar(int i)
{
/*
-----------------------------------------------------------------------

     spmpar provides the single precision machine constants for
     the computer being used.  It is assumed that the argument
     i is an integer having one of the values 1, 2, or 3.
     If the single precision arithmetic being used has m base b digits and
     its smallest and largest exponents are emin and emax, then

	spmpar(1) = b^(1 - m),		the machine precision,

	spmpar(2) = b^(emin - 1),	the smallest magnitude,

	spmpar(3) = b^emax*(1 - b^(-m)), the largest magnitude.

-----------------------------------------------------------------------
     written by
	Alfred H. Morris, Jr.
	Naval Surface Warfare Center (NSWC)
-----------------------------------------------------------------------
-----------------------------------------------------------------------
     Modified by BARRY W. BROWN to return DOUBLE PRECISION machine
     constants for the computer being used.  This modification was
     made as part of converting BRATIO to double precision
-----------------------------------------------------------------------
*/
    static double b,binv,bm1,one,w,z;
    static int emax,emin,ibeta,m;

    if(i == 1) {
	b = ipmpar(4);
	m = ipmpar(8);
	return(pow(b,(double)(1-m)));
    }
    if(i == 2) {
	b = ipmpar(4);
	emin = ipmpar(9);
	one = 1;
	binv = one/b;
	w = pow(b,(double)(emin+2));
	return(w*binv*binv*binv);
    }
    /* else  i = 3 : */
    ibeta = ipmpar(4);
    m = ipmpar(8);
    emax = ipmpar(10);
    b = ibeta;
    bm1 = ibeta-1;
    one = 1;
    z = pow(b,(double)(m-1));
    w = ((z-one)*b+bm1)/(b*z);
    z = pow(b,(double)(emax-2));
    return(w*z*b*b);
}

double stvaln(double *p)
{
/*
**********************************************************************

		    STarting VALue for Neton-Raphon
		calculation of Normal distribution Inverse

			      Function


     Returns X	such that CUMNOR(X)  =	 P,  i.e., the	integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


			      Arguments


     P --> The probability whose normal deviate is sought.
		    P is DOUBLE PRECISION


			      Method

     The  rational   function	on  page 95    of Kennedy  and	Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980.

**********************************************************************
*/
static double xden[5] = {
    0.993484626060e-1,0.588581570495,0.531103462366,0.103537752850,
    0.38560700634e-2 };
static double xnum[5] = {
    -0.322232431088,-1.000000000000,-0.342242088547,-0.204231210245e-1,
    -0.453642210148e-4 };

    double sign,y,z;

    if(*p > 0.5) {	sign =	1;	z = 1-*p; }
    else {		sign = -1;	z = *p; }
    y = sqrt(-(2*log(z)));
    return sign*(y + devlpl(xnum,5,y) / devlpl(xden,5,y));
}

/*----------------------------------------------------------------------
FIFDINT:
Truncates a double precision number to an integer and returns the
value in a double.
----------------------------------------------------------------------**/
double fifdint(double a)
/* a	 -     number to be truncated */
{
  long temp;
  temp = (long)(a);
  return (double)(temp);
}

/*----------------------------------------------------------------------**
FIFDSIGN:
transfers the sign of the variable "sign" to the variable "mag"
----------------------------------------------------------------------**/
double fifdsign(double mag,double sign)
/* mag	   -	 magnitude */
/* sign	   -	 sign to be transfered */
{
  if (mag < 0) mag = -mag;
  if (sign < 0) mag = -mag;
  return mag;

}
/*----------------------------------------------------------------------**
FIFIDINT:
Truncates a double precision number to a long integer
----------------------------------------------------------------------**/
long fifidint(double a)
/* a - number to be truncated */
{
  return (long)(a);
}

/*----------------------------------------------------------------------**
FTNSTOP:
Prints msg to standard error and then exits
----------------------------------------------------------------------**/
void ftnstop(char* msg)
/* msg - error message */
{
  if (msg != NULL) fprintf(stderr,"%s\n",msg);
  exit(EXIT_FAILURE); /* EXIT_FAILURE from stdlib.h, or use an int */
}
