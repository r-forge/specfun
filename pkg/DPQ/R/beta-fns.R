#### -- $Id: beta-fns.R,v 1.17 2014/05/28 14:55:58 maechler Exp maechler $

##- if(exists('lbeta', mode='function', inherits=F)) rm(lbeta)
##- lbeta0 <- lbeta

## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")

##' Qab := Q(a,b) := Gamma(a+b) / Gamma(b)  =  Gamma(a) / Beta(a,b)
##' <==>
##' log(Q(a,b)) = logQab(a,b)
##'             = lgamma(a+b) - lgamma(b) =
##'             = lgamma(a)   - lbeta(a,b)
logQab.asy <- function(a,b, k.max = 5, give.all = FALSE)
{
  ## Purpose:  log(Q(a,b)) asymptotically when max(a,b) -> Inf  &  a^2/b -> C
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  8 Jun 97, 17:28

  if(any(a<=0) || any(b<=0)) stop("a & b must be > 0")
  if(all(a >= b)) { cc <- b; b <- a; a <- cc } #- Now: a <= b
  else if(any(a>b)) stop('Require a <= b  or  a > b   for  ALL vector elements')
  if(length(b)>1) stop("this function only works for SCALAR b (> a)")
  na <- length(a)
  pa <- rbind(Qab.terms(a, k.max)) ##  --> pa = matrix   na  x  (k.max+1)
  ##          ~~~~~~~~~
  ki <- if(k.max >= 1) k.max:1 else numeric(0)
  if(give.all) {
    r <- matrix(0, nrow=na, ncol=k.max+1,
                dimnames = list(format(a), paste(k.max:0)))
    for(k in ki) r[,k] <- pa[,k] + r[,k+1]*(a-k-1)/b
  } else {
    r <- numeric(na)
    for(k in ki) r <- pa[,k] + r*(a-k-1)/b
  }
  a*log(b) + log(1 + a*(a-1)/b * r)
}

Qab.terms <- function(a, k)
{
  ## Purpose: Compute the terms used for  Qab, and beta function
  ## 		Qab := Q(a,b) := Gamma(a+b) / Gamma(b)
  ## -------------------------------------------------------------------------
  ## Arguments: a: the smaller of the arguments of  beta(a,b)
  ##            k: the number of terms in the series expansion
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 19 Jun 97, 11:08

  k <- as.integer(k)
  if(!is.numeric(a)) a <- as.numeric(a) #- leave integers alone..
  if(length(k)!=1 || k > 5 || k < 0)
    stop("'k' must be 1 number in { 0, 1,..,5 }")
  na <- length(a)
  Qab.trms <- expression(
      rep.int(1/2, na),	##  k=1
      a/8 - 1/24, 	##  k=2
      a*(a-1)/48,	##  k=3
      (2/5 + a*(1+3*a*(a-2)))/1152,  ## k=4
      a*(1 + a*(2+a*(a-5)))/5760     ## k=5
  )[seq_len(k)]

  vapply(Qab.trms, eval, numeric(na), envir=environment())
}

### NOTA BENE: All these  lbeta*() are -- I think -- not better than the
##  ---------  builtin  lbeta()  [which exists only since 1998]

lbeta.asy <- function(a,b, k.max = 5, give.all = FALSE)
{
  ## Purpose:  log(beta( a, b))  which also works for HUGE  a / b or b/a
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ##  log(B(a,b)) = log(G(a) G(b) / G(a+b)) =
  ##              = log(G(a)) - log( G(a+b)/G(b) ) = log(G(a)) - log(Q(a,b))
  lgamma(a) - logQab.asy(a,b, give.all = give.all)
}

## MM: this is unused (why ??)
lbetaMM <- function(a,b, Cut.asy = 1e-2)
{
  ## Purpose:  log(beta( a, b))  which also works for HUGE  a / b or b/a
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 May 97, 17:17
  if(a <= 0 || b <= 0) stop("a & b must be > 0")
  if(length(a)>1 || length(b) > 1)
    stop("this function only works for SCALAR arguments")
  stopifnot(is.numeric(Cut.asy), length(Cut.asy) == 1, Cut.asy >= 0)
  if(a>b) { cc <- b; b <- a; a <- cc } #- Now: a <= b
  if(a*a < b*Cut.asy) {
    if(exists("DEBUG") && DEBUG)
      warning(paste("a=",formatC(a)," b=",formatC(b),
                    " -- using asymptotic  lbeta(.)"))
    lgamma(a) - logQab.asy(a,b)
  } else
    lgamma(a)+lgamma(b)-lgamma(a+b) ## was =: lbeta00(a, b)
}

lbetaM <-  function(a,b, k.max = 5, give.all = FALSE)
{
  ## Purpose:  log(beta( a, b))  which also works for HUGE  a / b or b/a
  ## -------------------------------------------------------------------------
  ## Arguments: a,b: as for  lbeta(.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 May 97, 17:17
  if(any(a<=0) || any(b<=0)) stop("a & b must be > 0")
  na <- length(a)
  nb <- length(b)
  if(nb == 1 || nb == na) {
    r <- numeric(na)
    for(i in 1:na) r[i] <- lbeta.asy(a[i], b[if(nb>1)i else 1],
                                     k.max= k.max, give.all = give.all)
  } else { ##--- return  outer(..)
    r <- matrix(0, nrow = nb, ncol = na)
    for(i in 1:na)r[,i] <- lbeta.asy(a[i], b, k.max= k.max, give.all = give.all)
  }
  r
}
lbeta.n <- function(n, a)
{
  ## Purpose:  log(beta(a, n)) for (small) INTEGER n.
  ##    beta(a,n) = G(a) G(n) / G(a+n) = (n-1)! / (a*(a+1)*...*(a+n-1))
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 May 97, 10:08
  if(is.na(n <- as.integer(n)) || min(n) <= 0 || max(n) > 10000)
    stop("'n' must be positive (not too large) integer")
  r <- numeric(nn <- length(n))
  for (i in 1:nn)
    r[i] <- lgamma(n[i]) - sum(log(a+0:(n[i]-1)))
  r
}

beta.n<- function(n, a)
{
  ## Purpose: beta(a, n) for (small) INTEGER n.
  ##    beta(a,n) = G(a) G(n) / G(a+n) = 1*2*...*(n-1) / (a*(a+1)*...*(a+n-1))
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 May 97, 10:08
  if(length(n) > 1) stop("'n' must have length 1.")
  if(is.na(n <- as.integer(n)) || n <= 0)
    stop("'n' must be positive (not too large) integer")
  ni <- seq(length= n-1)
  prod(ni) / prod(a+c(0,ni))
}

G.half <- 1.7724538509055160272981674833411451827974 # sqrt(pi) = Gamma(1/2)

lbeta.n.half <- function(n, a)
{
  ## Purpose:  log(beta(a, n + 1/2)) for (small) INTEGER n.
  ##    beta(a,n + h) = G(a) G(n+h) / G(a+n+h) = (n-1)! / (a*(a+1)*...*(a+n-1))
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 May 97, 10:08
  if(length(n) > 1) stop("'n' must have length 1.")
  if(is.na(n <- as.integer(n)) || n <= 0)
    stop("'n' must be positive (not too large) integer")
  #............
  stop("unfinished -- FIXME")
}

### pbeta()  --- is now used from TOM708 exclusively
### -------  ~/R/D/r-devel/R/src/nmath/toms708.c
###                                    ~~~~~~~~~

## From ~/R/D/r-devel/R/src/nmath/pgamma.c :
##                                --------
##/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */

SQR <- function(x) (x*x)
scalefactor <- SQR(SQR(SQR(4294967296.0)))
rm(SQR)

##/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ] =~=  -x */
M_cutoff <- log(2) *
    .Machine$double.max.exp / .Machine$double.eps ## = 3.196577e18

logcf <- function (x, i, d, eps, maxit=10000) ##/* ~ relative tolerance */)
{
    ## Continued fraction for calculation of  sum_{k=0}^Inf x^k/(i+k*d) =
    ##		1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ...
    ##
    ## auxiliary in log1pmx() and lgamma1p()

    stopifnot (i > 0, d >= 0, eps > 0)

    c1 <- 2 * d
    c2 <- i + d
    c4 <- c2 + d
    a1 <- c2
    b1 <- i * (c2 - i * x)
    ddx <- d * d * x
    A2 <- c4 * c2 - ddx
    B2 <- c4 * b1 - i * ddx

    it <- 1## vectorize in x (!)
    while (any(abs(A2 * b1 - a1 * B2) > abs(eps * b1 * B2))) {
        c3 <- c2*c2*x
	c2 <- c2 + d
	c4 <- c4 + d
	a1 <- c4 * A2 - c3 * a1
	b1 <- c4 * B2 - c3 * b1

	c3 <- c1 * c1 * x
	c1 <- c1 + d
	c4 <- c4 + d
	A2 <- c4 * a1 - c3 * A2
	B2 <- c4 * b1 - c3 * B2

        m.B2 <- exp(mean(log(abs(B2))))
	if (m.B2 > scalefactor) {
	    a1 <- a1 / scalefactor
	    b1 <- b1 / scalefactor
	    A2 <- A2 / scalefactor
	    B2 <- B2 / scalefactor
	} else if (m.B2 < 1 / scalefactor) {
	    a1 <- a1 * scalefactor
	    b1 <- b1 * scalefactor
	    A2 <- A2 * scalefactor
	    B2 <- B2 * scalefactor
	}
        it <- it+1
        if(it > maxit)
            warning("non-convergence in logcf(), iter > ", maxit)
    }
    ## return
    (A2 / B2)
}


## Accurate calculation of log(1+x)-x, particularly for small x.
log1pmx <- function(x, tol_logcf = 1e-14) {
    minLog1Value <- -0.79149064

    r <- x
    if(any(c1 <- (x > 1 | x < minLog1Value)))
        r[c1] <- log1p(x[c1]) - x[c1]
    ## else { ## ##/* expand in [x/(2+x)]^2 */
    if(any(c2 <- !c1)) {
        x <- x[c2]
	term <- x / (2 + x)
	y <- term * term
        r[c2] <- term *
            ifelse(abs(x) < 1e-2,
                   (((2 / 9 * y + 2 / 7) * y + 2 / 5) * y + 2 / 3) * y - x,
                   2 * y * logcf(y, 3, 2, tol_logcf) - x)
    }
    r
}

lgamma1p <- function(a, tol_logcf = 1e-14)
{
    ## Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5).

    eulers_const <- 0.5772156649015328606065120900824024

    ##/* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    N <- 40
    coeffs <-
        c(0.3224670334241132182362075833230126e-0, ##/* = (zeta(2)-1)/2 */
          0.6735230105319809513324605383715000e-1, ##/* = (zeta(3)-1)/3 */
          0.2058080842778454787900092413529198e-1,
          0.7385551028673985266273097291406834e-2,
          0.2890510330741523285752988298486755e-2,
          0.1192753911703260977113935692828109e-2,
          0.5096695247430424223356548135815582e-3,
          0.2231547584535793797614188036013401e-3,
          0.9945751278180853371459589003190170e-4,
          0.4492623673813314170020750240635786e-4,
          0.2050721277567069155316650397830591e-4,
          0.9439488275268395903987425104415055e-5,
          0.4374866789907487804181793223952411e-5,
          0.2039215753801366236781900709670839e-5,
          0.9551412130407419832857179772951265e-6,
          0.4492469198764566043294290331193655e-6,
          0.2120718480555466586923135901077628e-6,
          0.1004322482396809960872083050053344e-6,
          0.4769810169363980565760193417246730e-7,
          0.2271109460894316491031998116062124e-7,
          0.1083865921489695409107491757968159e-7,
          0.5183475041970046655121248647057669e-8,
          0.2483674543802478317185008663991718e-8,
          0.1192140140586091207442548202774640e-8,
          0.5731367241678862013330194857961011e-9,
          0.2759522885124233145178149692816341e-9,
          0.1330476437424448948149715720858008e-9,
          0.6422964563838100022082448087644648e-10,
          0.3104424774732227276239215783404066e-10,
          0.1502138408075414217093301048780668e-10,
          0.7275974480239079662504549924814047e-11,
          0.3527742476575915083615072228655483e-11,
          0.1711991790559617908601084114443031e-11,
          0.8315385841420284819798357793954418e-12,
          0.4042200525289440065536008957032895e-12,
          0.1966475631096616490411045679010286e-12,
          0.9573630387838555763782200936508615e-13,
          0.4664076026428374224576492565974577e-13,
          0.2273736960065972320633279596737272e-13,
          0.1109139947083452201658320007192334e-13 ##/* = (zeta(40+1)-1)/(40+1) */
          )

    c <- 0.2273736845824652515226821577978691e-12##/* zeta(N+2)-1 */

    r <- a
    ## if(abs(a) >= 0.5)
    a.lrg <- (abs(a) >= 0.5)
    r[a.lrg] <- lgamma(a[a.lrg] + 1)
    ## else {
    if(any(a.sml <- !a.lrg)) {
        ##/* Abramowitz & Stegun 6.1.33,
        ## * also  http://functions.wolfram.com/06.11.06.0008.01 */
        a. <- a[a.sml]
        lgam <- c * logcf(-a. / 2, N + 2, 1, tol_logcf)
        for (i in N:1)
            lgam <- coeffs[i] - a. * lgam

        ## return
        r[a.sml] <- (a. * lgam - eulers_const) * a. - log1pmx(a.)
    }
    r
}## lgamma1p




##' original (getting inaccurate when p --> 1)
qnorm.appr <- function(p) {
  ## qnorm: normal quantile approximation for qnorm(p),  for p > 1/2
  ## -- to be used in  qbeta(.)
  ## The relative error of this approximation is quite ASYMMETRIC: mainly < 0
  C1 <- 2.30753
  C2 <- 0.27061
  C3 <- 0.99229
  C4 <- 0.04481
  r <- sqrt(-2*log(1-p))
  r - (C1 + C2 * r) / (1 + (C3 + C4 * r) * r)
}

qnormU.appr <- function(p,
                        lp = R.DT.Clog(p, lower.tail=lower.tail, log.p=log.p),
                                        # ~= log(1-p) -- independent of lower.tail, log.p
                        lower.tail=TRUE, log.p=FALSE)
{
    ## qnorm: normal quantile approximation;
    ## -- to be used in  qbeta(.)
    ## The relative error of this approximation is quite ASYMMETRIC: mainly < 0
    C1 <- 2.30753
    C2 <- 0.27061
    C3 <- 0.99229
    C4 <- 0.04481

    if(missing(p)) { ## have lp = log(1-p)  <==>  e^lp = 1-p  <==>  p = 1 - e^lp
        p. <- -expm1(lp)
        ## swap p <--> 1-p -- so we are where approximation is better
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector

    } else {
        p. <- R.D.qIv(p, log.p)
        ## swap p <--> 1-p -- so we are where approximation is better
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector
        p[swap] <- if(log.p) R.Log1.Exp(p[swap]) else 1 - p[swap]
    }
    R <- r <- sqrt(-2 * lp)
    r.ok <- r < 1e10 ## e.g., for p == 1, r = Inf would give R = NaN
    ## for r >= 1e10:  R = r numerically in formula below
    r <- r[r.ok]
    R[r.ok] <- r - (C1 + C2 * r) / (1 + (C3 + C4 * r) * r)
    R[swap] <- -R[swap]
    R[p. == 1/2] <- 0
    R
}

qbeta.appr.1 <- function(a, p, q, y = qnormU.appr(a, lower=FALSE))
{
  ## Purpose: Approximate  qbeta(a, p,q) -- Abramowitz & Stegun (26.5.22)
  ##          qbeta(.) takes this only when  p>1 & q>1
  ## -------------------------------------------------------------------------
  ## Arguments: a: percentage point;  (p,q): beta parameters
  r <- (y * y - 3) / 6 ## lambda
  s <- 1 / (p + p - 1)
  t <- 1 / (q + q - 1)
  h <- 2 / (s + t)
  w <- y * sqrt(h + r) / h - (t - s) * (r + 5 / 6 - 2 / (3 * h))
  p / (p + q * exp(w + w))
}

qbeta.appr.3 <- function(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
{
    ##  a=alpha
    ## Purpose: Approximate  qbeta(a, p,q) -- for small a --
    ##  Inversion of   I_x(a,b) ~= x^a / (a B(a,b))

    ## CARE:  log(p) + logbeta == log(p * beta(p,q)) suffers from cancellation
    ## ----   for small p
    ## ---> look at Qab() above or *rather* also see experiments c12pBeta() , p.err.pBeta()
    ##  in  ./qbeta-dist.R  <<<<
    ##        ============  <<<<

    ## log
    log.a <- if(lower.tail) R.D.log(a, log.p=log.p) else R.D.LExp(a, log.p=log.p)

    pmin(1, exp(log.a +  log(p) + logbeta) / p)
}

qbeta.appr.2 <- function(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
{
    ## Purpose: Approximate  qbeta(a, p,q)

    l1ma <- if(lower.tail) R.D.LExp(a, log.p=log.p) else R.D.log(a, log.p=log.p)
    ## pmax(0, ....) -- needed when  (l1ma + log(q) + logbeta) / q > 0, e.g., for
    ## e.g. qbeta.appr.2(1/4, 1/2, 1/2)
    -expm1((l1ma + log(q) + logbeta) / q)
}

qbeta.appr.4 <- function(a, p, q, y = qnormU.appr(a, lower=FALSE),
                         verbose=getOption("verbose"))
{
    ## Purpose: Approximate  qbeta(a, p,q); 'a' is only used via qnorm(a,..)
    r <- q + q
    t <- 1 / (9 * q)
    t <- r * (1 - t + y * sqrt(t))^3 # = \chi^2_alpha .. which "must be > 0", but I find: no!
    if(verbose) cat("chi^2[alpha] =: t = ",format(t), if(t <= 0)"t < 0 !!","\n")
    t <- (4 * p + r - 2) / t
    if(verbose) cat("t = (1 + x0)/(1 - x0) = ",format(t), if(t <= 1) "t <= 1 !!","\n")
    1 - 2 / (t + 1)
}

qbeta.appr <- function(a, p, q, y = qnormU.appr(a, lower=FALSE), logbeta= lbeta(p,q),
                       verbose = getOption("verbose") && length(a) == 1)
{
    ## Purpose: Approximate  qbeta(a, p,q) --- for  a <= 1/2
    ## -------------------------------------------------------------------------
    ## Arguments:
    ## -------------------------------------------------------------------------
    ## Author: After qbeta.c in  R  (--> AS 109...)
    if(length(p)!=1 || length(q)!=1)
        stop("'p' & 'q' must have length 1 !!")
    if(any(a > 1/2)) warning("a[.] > 1/2 (not thought for!)")

    if (p > 1 && q > 1) { ## Abramowitz & Stegun(26.5.22)
        if(verbose) cat("p,q > 1: Using qbeta.appr.1(): Abramowitz-Stegun (26.5.22)\n")
        ## 'a' is not needed here, just 'y' is
        qbeta.appr.1(,p,q,y)
    } else {
        if(verbose) cat("p or q <= 1 : ")
        r <- q + q
        ## more "robustly"
        ## t <- 1 / (9 * q)
        ## t <- r * (1 - t + y * sqrt(t))^3
        st <- 1/(3 * sqrt(q))
        t <- r * (1 - st*(st + y))^3

        ## NOTE: length(t) == length(y) == length(a)  is fulfilled

        ## vectorized in t { ==> logical vectors instead of simple if(..) .. else if(..) .. else ..
        ans <- t
        if(any(neg <- t <= 0)) { ##- forget t(q, y)  and  y(a)
            if(verbose)   cat("t <= 0:         appr. 2 (beta, q)\n")
            ans[neg] <- - expm1((log1p(-a[neg]) + log(q) + logbeta) / q)
        }
        if(any(!neg)) { ##else
            t <- (4 * p + r - 2) / t
            if (any(L1 <- !neg & t <= 1)) {
                if(verbose) cat("t > 0; t' <= 1: appr. 3 (beta, p)\n")
                ## ans[L1] <- exp((log(a[L1] * p) + logbeta) / p)
                ans[L1] <- qbeta.appr.3(a[L1], p, q, logbeta = logbeta)
                ## ==> linear relationship on log-log scale:
                ##   log(ans) = 1/p * log(a) + (log p + log beta(p,q))/p
            }
            if (any(L2 <- !neg & t > 1)) { # the other case
                if(verbose) cat("t > 0; t' > 1 : appr. 4 = 1-2/(t'+1)\n")
                ans[L2] <- 1 - 2 / (t[L2] + 1)
            }
        }
        ans
    }
}

## Hmm, this is (sometimes!) not quite equivalent to C's :
## ex.  qbeta.R(0.5078, .01, 5) -> 2.77558e-15    but qbeta() gives 1.776357e-15
qbeta.R	 <-  function(alpha, p, q,
                      lower.tail = TRUE, log.p = FALSE,
		      logbeta = lbeta(p,q),
		      low.bnd = 3e-308, up.bnd = 1-2.22e-16,
                      method = c("AS109", "Newton-log"),
                      ## FIXME: (a,p,q) : and then uses (a, pp) .. hmm
		      f.acu = function(a,p,q) max(1e-300, 10^(-13- 2.5/pp^2 - .5/a^2)),
		      fpu = .Machine$ double.xmin,
		      qnormU.fun = function(u, lu) qnormU.appr(p=u, lp=lu, lower=FALSE),
                      R.pre.2014 = FALSE,
		      verbose = getOption("verbose") ## FALSE, TRUE, or 0, 1, 2, ..
		      )
{
### ----- following the C code in  ~/R/D/r-devel/R/src/nmath/qbeta.c
###
### ----------> make use of 'log.p' in	qnormU.appr(* , log.p)	as well!

    ## Purpose:	 qbeta(.) in R -- edited;  verbose output
    ## -------------------------------------------------------------------------
    ## Arguments: alpha: percentage;  (p,q) Beta parameter
    ##	  low.bnd, up.bnd: algorithm parameter
    ##	acu, fpu:	accuracy;     .Machine$ double.xmin = 2.22e-308 for IEEE
    ##		SAE below is the most negative decimal exponent which does not
    ##		cause an underflow; a value of -308 or thereabouts will often be
    ##		OK in double precision.
    ##		sae <- -37 ;  fpu <- 10 ^ sae
    ##	  qnormU.fun(p) ~= qnorm(1 - p) = qnorm(p, ..., lower=FALSE)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Apr 1997, 12:55
    ## ------ edited C-code  qbeta.c  into  R-code

    if(length(alpha)>1) stop("only works for length(1) 'alpha'")

    if (low.bnd >= up.bnd || low.bnd < 0 || up.bnd > 1)
	stop("MUST  0 <= low.bnd < up.bnd <= 1")

    ## Test for admissibility of parameters
    if (p < 0 || q < 0 || alpha < 0.0 || alpha > 1.0) stop("DOMAIN_ERROR")

    ##  R_Q_P01_boundaries(alpha, 0, 1) :=
    if(log.p) {
        if(alpha == -Inf) return(if(lower.tail) 1 else 0)
        if(alpha ==   0 ) return(if(lower.tail) 0 else 1)
    } else {
        if (alpha == 0 || alpha == 1) return(if(lower.tail) alpha else 1-alpha)
    }

    p. <- R.DT.qIv(alpha, lower.tail, log.p=log.p) # = lower_tail prob (in any case)

    if(R.pre.2014) {
        if(log.p && (p. == 0. || p. == 1.))
	return(p.) ## better than NaN or infinite loop;

    } else { ## can and *must* do better:
        if(log.p && (p. == 0. || p. == 1.))
            warning("p. is 0 or 1")

        if(FALSE) {
            lp <- -10*(150:20); qb <- qbeta(lp, 2,3, log=TRUE)
            plot(lp, qb, ylab = "qbeta(lp, 2,3, log=TRUE)", log="y", type="o")
            cbind(lp, qb, pbeta(qb, 2,3, log=TRUE))
        }
    }

    ## inflection point of pbeta() (outside [0,1] when p < 1 < q or q < 1 < p)
    x.ip <- (p-1) / (p+q -2)
    y.ip <- pbeta(x.ip, p,q)

    ## change tail if necessary, such that  beta(a,p,q)	 with  0 = a <= 1/2
    ## la := log(a), but without numerical cancellation:
    if (alpha <= 1/2) {
	a <- p.
	la <- if(lower.tail) R.D.log(alpha, log.p) else R.D.LExp(alpha, log.p)
        pp <- p; qq <- q; swap.tail <- FALSE
    } else {
	a <- R.DT.CIv(alpha, lower.tail, log.p)
	la <- if(lower.tail) R.D.LExp(alpha, log.p) else R.D.log(alpha, log.p)
        pp <- q; qq <- p; swap.tail <- TRUE
    }
    if(verbose) cat(if(swap.tail)"SWAP tail" else "no swap", "\n")
    stopifnot( all.equal(log(a), la) )

    acu <- f.acu(a,p,q)
    ## set the exponent of acu to -2r-2 for r digits of accuracy
    if(acu <= 0 || acu > .1) acu   <- 1.0e-32 # 32: 15 digits accu.

    if(verbose) cat("logbeta=", formatC(logbeta,dig=19),
		    "; acu=",formatC(acu, dig=6) ,"\n")

    ## calculate the initial approximation 'xinbta' [FIXME: never designed for log.p=TRUE]
    x0 <- xinbta <- qbeta.appr(a, pp,qq, y = qnormU.fun(u=a, lu=la), logbeta=logbeta,
                               verbose = verbose)

    if(verbose) cat("initial 'xinbta' =", formatC(xinbta, dig=18))

    xinbta0 <- xinbta
    if(!is.finite(xinbta)) {
	warning("** qbeta.R(): qbeta.appr() was NOT finite!! **")
        xinbta <- 0.5
    }

    ## Solve for x by a modified Newton-Raphson method,
    ## using the function  pbeta(.)

    r <- 1 - pp
    t <- 1 - qq
    yprev <- 0
    adj <- 1
    ##not used: prev <- 1
    if (xinbta < low.bnd) { xinbta <- 0.5; if(verbose)cat(" -- low.bnd  -> xi= 0.5\n") } else
    if (xinbta >  up.bnd) { xinbta <- 0.5; if(verbose)cat(" -- up.bnd   -> xi= 0.5\n") } else
    if(verbose) cat("\n")

    ##-- for single precision -- from CORRECTED	 AS 109:
    ##R	 iex <- as.integer( max(-5/pp^2 - 1/a^2 - 13, sae))
    ##NN acu <- max(fpu, 10 ^ (-5/pp^2 - 1/a^2 - 13))

    finish <- FALSE
    o.it <- 0
    method <- match.arg(method)

    if(method == "AS109")
    while(!finish) { ##-- Outer Loop
	y <- pbeta(xinbta, pp, qq)
	if(verbose) cat("y=pbeta(x..)=", formatC(y, dig=15,wid=18, flag='-'))
	if(!is.finite(y)) {
	    cat("  y=pbeta(): not finite: 'DOMAIN_ERROR'\n",
		"  xinbta =",form01.prec(xinbta,dig=18),"\n")
	    return(xinbta) ##browser()
	}
	y <- (y - a) * exp(logbeta + r * log(xinbta) + t * log1p(-xinbta))
	if(verbose) cat(" y.n=(y-a)*e^..=", formatC(y))
	if(!is.finite(y)) {
	    cat("  y = (y-a)*exp(...) not finite:\n",
		"  xinbta =",form01.prec(xinbta,dig=18),"\n")
 	    return(xinbta) ##browser()
	}

	if (y * yprev <= 0) prev <- max(fpu, abs(adj))
	g <- 1
	in.it <- 0
	if(verbose>=2) cat("\n	Inner loop: ")
	for(in.it in 1:1000) {
	    adj <- g * y
	    if (abs(adj) < prev) { #-- current adjustment < last one which gave sign change
		if(verbose>=2)cat("<")
		tx <- xinbta - adj
		if (0 <= tx && tx <= 1) {
		    if (prev <= acu || abs(y) <= acu) { ## goto L_converged
			finish <- TRUE; break
		    }
		    if (tx != 0 && tx != 1) break
		}
	    } else if(verbose>=2) cat(".")
	    if(!finish) g <- g / 3
	} ##-- end inner loop
	if(verbose>=2) cat(if(in.it>10)"\n" else " ")
	if(verbose) cat(" ",in.it,"it --> adj=g*y=",formatC(adj),"\n")
	if (abs(tx - xinbta) < 1e-15*xinbta) break # goto L_converged;
	if(verbose>=2) cat("--Outer Loop: |tx - xinbta| > eps; tx=",formatC(tx, dig=18),
	   " tx-xinbta =", formatC(tx-xinbta),"\n")
	xinbta <- tx
	yprev <- y
	o.it <- o.it + 1
    } ##- end outer loop
    else if(method == "Newton-log") #--------------------------------- new --------
        ## but really, I want  Newton on   log-*log* scale --
    while(!finish) { ##-- Outer Loop
	y <- pbeta(xinbta, pp, qq, log = TRUE)
	if(verbose) cat("y=pbeta(xi,*,log=TRUE)=", formatC(y, dig=15,wid=18, flag='-'))
	if(!is.finite(y)) {
	    cat("  y=pbeta(): not finite: 'DOMAIN_ERROR'\n",
		"  xinbta =",formatC(xinbta,dig=18),"\n")
	    return(xinbta) ##browser()
	}
        ## y := g(xinbta) / g'(xinbta)  where g(x) =  log(alpha) - log pbeta(x, ..)
        ##                              ==>  g'(x) = - dbeta(x,.) / pbeta(x,.)
	y <- (y - la) * exp(y + (logbeta + r * log(xinbta) + t * log1p(-xinbta)))
	if(verbose) cat(" y.n=(y-a)*e^..=", formatC(y))
	if(!is.finite(y)) {
	    cat("  y = (y-a)*exp(...) not finite:\n",
		"  xinbta =",formatC(xinbta,dig=18),"\n")
 	    return(xinbta) ##browser()
	}

## 	if (y * yprev <= 0) prev <- max(fpu, abs(adj))
## 	g <- 1
## 	in.it <- 0
## 	if(verbose>=2) cat("\n	Inner loop: ")
## 	for(in.it in 1:1000) {
## 	    adj <- g * y
## 	    if (abs(adj) < prev) { #-- current adjustment < last one which gave sign change
## 		if(verbose>=2)cat("<")
## 		tx <- xinbta - adj
## 		if (0 <= tx && tx <= 1) {
## 		    if (prev <= acu || abs(y) <= acu) { ## goto L_converged
## 			finish <- TRUE; break
## 		    }
## 		    if (tx != 0 && tx != 1) break
## 		}
## 	    } else if(verbose>=2) cat(".")
## 	    if(!finish) g <- g / 3
## 	} ##-- end inner loop
## 	if(verbose>=2) cat(if(in.it>10)"\n" else " ")
## 	if(verbose) cat(" ",in.it,"it --> adj=g*y=",formatC(adj),"\n")
 	if(verbose) cat("\n")
        tx <- xinbta - y
	if (abs(tx - xinbta) < 1e-15*xinbta) break # goto L_converged;
	if(verbose>=2) cat("--Outer Loop: |tx - xinbta| > eps; tx=",formatC(tx, dig=18),
	   " tx-xinbta =", formatC(tx-xinbta),"\n")
	xinbta <- tx
## 	yprev <- y
	o.it <- o.it + 1
    } ##- end outer loop
    else stop("unknown method ", dQuote(method))

    ## L_converged:
    if(verbose) cat(" ", o.it,"outer iterations\n")
    if (swap.tail) 1 - xinbta else xinbta
}

### form01.prec  -> ../../../MISC/Util.R
## source("~/R/MM/MISC/Util.R")
form01.prec <- function(x, digits = getOption("digits"), width = digits + 2,
                        eps = 1e-6, ...,
                        fun = function(x,...) formatC(x, flag='-',...))
{
  ## Purpose: format numbers in [0,1] with "precise" result,
  ##          using "1-.." if necessary.
  ## -------------------------------------------------------------------------
  ## Arguments: x:      numbers in [0,1]; (still works if not)
  ##            digits, width: number of digits, width to use with 'fun'
  ##            eps:    Use '1-' iff  x in  (1-eps, 1] -- 1e-6 is OPTIMAL
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14 May 97, 18:07
  if(as.integer(digits) < 4) stop('digits must be >= 4')
  if(eps < 0 || eps > .1) stop('eps must be in [0, .1]')
  i.swap <- 1-eps < x  &  x <= 1 #-- Use "1- ." if <~ 1,  normal 'fun' otherwise
  r <- character(length(x))
  if(any(i.swap))
    r[i.swap] <-
      paste("1-", fun(1-x[i.swap], digits=digits - 5, width=width-2, ...),
            sep='')# -5: '1-' + 4 for exponent -1 for '0' (in other case)
  if(any(!i.swap))
    r[!i.swap] <- fun(x[!i.swap], digits=digits, width=width,...)
  attributes(r) <- attributes(x)
  r
}

