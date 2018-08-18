#### R functions for the -- noncentral chisquare density -- dchisq(.,., ncp > 0)

### ---> ./chisq-nonc-ex.R
###      ~~~~~~~~~~~~~~~~~

## source("~/R/MM/NUMERICS/dpq-functions/dpq-h.R")# the macros

## R's builtin  dchisq(x, df, ncp) is now {since P.Dalgaards improvement}
## -----------                     definitely better than this:
dnoncentchisq <- function(x, df, del)
{
    ## The noncentral-chisquare density at x
    ## for df degrees of freedom, noncentrality parameter del.
    ## (Only) x may be a vector.
    ##
    if(length(del) > 1 || del < 0)
        stop("noncentrality parameter del must be scalar >= 0!")
    if(length(df) > 1 || df <= 0)
        stop("df must be scalar > 0!")
    upper <- floor(del/2 + 5 * (del/2)^0.5)
    kv <- 0:upper
    poiv <- if(del == 0) 1 else dpois(kv, del/2)
    n <- length(x)
    fv <- matrix(rep(0, n * (1 + upper)), n)
    for(k in kv)
        fv[, 1 + k] <- dchisq(x, 2 * k + df)
    c(fv %*% poiv)
}

dchisqAsym <- function(x, df, ncp, log = FALSE)
{
  ## Purpose: Asymptotic approximation of NON-central by central Chi^2
  ## ----------------------------------------------------------------------
  ## Arguments: as for dchisq()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  1 Apr 2008, 10:34
    nl <- df + ncp
    n2l <- nl + ncp
    ic <- nl/n2l # = 1/c = 1 / (1+b) -- 'b' of Abramowitz & Stegun, p.942 [26.4.27]
    ff <- dchisq(x*ic, df = nl*ic, log = log)
    if(log) log(ic) + ff else ic * ff
}

dnchisqBessel <- function(x, df, ncp, log = FALSE)
{
    ## Purpose: Implement Fisher(1928) = Johnson,Kotz & Bala. (29.4) [p.436]
    ##          == Exact Formula for non-central dchisq(), using besselI()
    ## ----------------------------------------------------------------------
    ## Arguments: as for dchisq()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr 2008, 21:09

    if (any(is.na(x)) || any(is.na(df)) || any(is.na(ncp)))
	return(x + df + ncp)
    if (any(ncp <= 0) || any(df < 0))
	stop("must have ncp > 0 (here), df >= 0")
    if (!all(is.finite(df)) || !all(is.finite(ncp)))
	stop("some (ncp, df) values are not finite")

    ## Recycle:
    n <- max(length(x), length(df), length(ncp))
    if(n > 1) {
        if(length(x) < n) x <- rep(x, length = n)
        if(length(df) < n) df <- rep(df, length = n)
        if(length(ncp) < n) ncp <- rep(ncp, length = n)
    }

    nu. <- (df - 2)/2
    y <- sqrt(ncp*x)
    ## Iy := besselI(y, ., expon.scaled=TRUE) = exp(-y) * I_.(y)
    Iy <- besselI(y, nu = nu., expon.scaled = TRUE)
    ## NOTE: Bessel(y) for large y --> ../bessel-large-x.R
    ## NOTE 2: besselI(..) also needs a  'log = FALSE' argument <<< !
    if(log) {
        (y -(ncp+x)/2) + (nu./2) * log(x/ncp) + log(Iy/2)
    }
    else { ## not log
        exp(y -(ncp+x)/2) * (x/ncp)^ (nu./2) * Iy/2
    }
}

p.dnchiB <- function(df, ncp, log=FALSE, from=0, to = 2*ncp, p.log="", ...)
{
    ## Purpose: Comparison plot of dchisq() and its Bessel-approximation
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr 2008, 21:40
    x <- NULL # -Wall(codetools)
    curve(dnchisqBessel(x, df=df, ncp=ncp, log=log), from=from, to=to,
          n=2001, ylab = "dnchisq[Bessel](x,*)",
          main = deparse(sys.call()), log=p.log, ...)
    curve(dchisq       (x, df=df, ncp=ncp, log=log), col=2,   n=4001,
          lty=3, lwd=2, add=TRUE)
    legend("topright", c("dchisq()", "dnchisqBessel()"),
           col = 2:1, lwd=2:1, lty=c(3,1), inset = .02)
}

###--- Now, the R version of the C implementation
## -->  ~/R/D/r-devel/R/src/nmath/dnchisq.c

dnchisqR <- function(x, df, ncp, log = FALSE)
{
    stopifnot(length(x) == 1, length(df) == 1,
              length(ncp) == 1, length(log) == 1,
              is.numeric(x), is.numeric(df),
              is.numeric(ncp), is.logical(log))

    if(is.na(x) || is.na(df) || is.na(ncp))
	return(x + df + ncp)

    if (ncp < 0 || df <= 0 || !is.finite(df) || !is.finite(ncp))
	return(NaN)

    if(x < 0) return(.D_0(log))
    if(x == 0 && df < 2.)
	return(Inf)
##     if(ncp == 0)
## 	return(dchisq(x, df, log=log))

    eps <- 5e-15
    ncp2 <- 0.5 * ncp

    ##/* find max element of sum */
    imax  <- ceiling((-(2+df) +sqrt((2-df) * (2-df) + 4 * ncp * x))/4)
    if (imax < 0)
        imax <- 0
    if(is.finite(imax)) {
        dfmid  <- df + 2 * imax
        ## mid = dpois_raw(imax, ncp2, FALSE) * dchisq(x, dfmid, FALSE)
        mid <- dpois(imax, ncp2, log=FALSE) * dchisq(x, dfmid, log=FALSE)
    } else mid <- 0

    if(mid == 0) {
	##/* underflow to 0 -- maybe numerically correct; maybe can be more accurate,
        ## particularly when  give_log = TRUE */
        ##/* Use  central-chisq approximation formula when appropriate;
        ##    * ((FIXME: the optimal cutoff also depends on (x,df);  use always here? )) */
        if(log || ncp > 1000.) {
            nl <- df + ncp; ic <- nl/(nl + ncp) ##/* = "1/(1+b)" Abramowitz & St.*/
            return(dchisq(x*ic, nl*ic, log=log))
        } else return(.D_0(log))
    }

    sum <- mid

    ##/* errorbound := term * q / (1-q)  now subsumed in while() / if() below: */

    ##/* upper tail */
    term <- mid; df <- dfmid; i <- imax
    repeat {
	i <- i+1
	q <-  x * ncp2 / i / df;
	df <- df+2
	term <- term*q
	sum <- sum + term
        if(!(q >= 1 || term * q > (1-q)*eps)) break
    }
    ##/* lower tail */
    term <- mid; df <- dfmid; i <- imax
    while(i) {
	df <- df - 2
	q <- i * df / x / ncp2
	i <- i-1
	term <- term * q
	sum <- sum + term
	if (q < 1 && term * q <= (1-q)*eps) break
    }
    ## return
    .D_val(sum, log)
}
