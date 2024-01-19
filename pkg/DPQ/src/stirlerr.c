/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000-2024, The R Core Team
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
 *
 *  DESCRIPTION
 *
 *    Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to
 *        write lgamma!)
 *
 * Merge in to R:
 * Copyright (C) 2000, The R Core Team
 * R has lgammafn, and lgamma is not part of ISO C
 */

#include "DPQpkg.h"

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 *
 * NB: stirlerr(n/2) & stirlerr((n+1)/2) are called from dt(x,n) for all real n > 0 ;
 *     stirlerr(n/2) from gamma(n/2) when n is integer and 10 < |n|/2 <= 50
 *     stirlerr(x)   from dpois_raw(x, lam) for any  x > 0; --> ./dpois.c [dpois_raw() called by many]
 *     stirlerr(n), stirlerr(x), stirlerr(n-x) from binom_raw(x, n, ..) for all possible 0 < x < n
 */

/* NB 2 --- in DPQ's R code, we have much improved and versatile stirlerr() ==> ../R/dgamma.R
 * --       This is *really* only for reproducibility of R's C level (non-API!) stirlerr()
_________ CURRENTLY *NOT* directly called from R, only used in gammafn_ver() ==> ./gamma-variants.c

_________ TODO:  add  'int version'  or 'order' or   double* cutoffs ________________
 */
double dpq_stirlerr(double n, stirlerr_version_t version)
{
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
/* From S5 on: ---> ../Misc/stirlerr-trms.R  ==> gmp::BernoulliQ() etc */
#define  S5   0.0019175269175269175269175262 // 691/360360
#define  S6   0.0064102564102564102564102561 // 1/156
#define  S7   0.029550653594771241830065352  // 3617/122400
#define  S8   0.17964437236883057316493850   // 43867/244188
#define  S9   1.3924322169059011164274315    // 174611/125400
#define S10  13.402864044168391994478957     // 77683/5796

/*
  exact values for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };

#define nC 10
    const static double cutoffs[nC] = {
	5.4,
	7.5,
	8.5,
	10.625,
	12.125,
	20,
	26,  //  15
	60,  //  35
	200, //  80
	3300 // 500
    };

    double nn;

  switch(version) { // stirlerr_version_t  <--> ./DPQpkg.h
  case R_3:   // "R<=3"
  case R_3_1: // "R4..1" : using lgamma1p()
    if (n <= 15.0) {
	nn = n + n;
	if (nn == (int)nn) return sferr_halves[(int)nn];
	return ((version == R_3) ? lgamma(n + 1.) : lgamma1p(n))
	    - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI;
    }
    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);

  case R_4_4: // "R4.4_0"  --- see R code in ../R/dgamma.R  <<<<<<
      nn = n + n;
      if (n <= 15. && (nn == (int)nn)) return sferr_halves[(int)nn];
      // else:
      if (n <= cutoffs[0])
	  return lgamma1p(n) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI;
      // else n > cutoffs[0] == 5.4
      nn = n*n;
      int k = nC;
      if (n > cutoffs[--k]) return (S0 -S1/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1 -S2/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2 -S3/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3 -S4/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4 -S5/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5 -S6/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6 -S7/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7 -S8/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8 -S9/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      /* else if (n > cutoffs[--k])*/
			    return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-S10/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;

  default: error(_("switch() logic programming errors in '%s()'"), "dpq_stirlerr");
  } // switch()
}

SEXP R_dpq_stirlerr(SEXP x_, SEXP stirlerr_v) {
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    R_xlen_t i, n = XLENGTH(x_);
    SEXP ans = PROTECT(allocVector(REALSXP, n));
    double *x = REAL(x_), *r = REAL(ans);
    stirlerr_version_t version = (stirlerr_version_t) asInteger(stirlerr_v);

    for(i=0; i < n; i++) {
	r[i] = dpq_stirlerr(x[i], version);
    }
    UNPROTECT(2);
    return ans;
}

