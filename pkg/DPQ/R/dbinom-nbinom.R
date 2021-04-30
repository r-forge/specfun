#### R version of  dbinom_raw()          from ~/R/D/r-devel/R/src/nmath/dbinom.c
####        and dnbinom() / dnbinom_mu()  "   ~/R/D/r-devel/R/src/nmath/dnbinom.c
##
##  (orig. was ~/R/MM/NUMERICS/dpq-functions/dbinom_raw.R )

dbinom_raw <- function(x, n, p, q=1-p, log = FALSE # >> ../man/dbinom_raw.Rd <<<
                       , verbose = getOption("verbose")
                       )
{
  ## Purpose: R version of dbinom_raw()   { from .../R/src/nmath/dbinom.c }
  ## ----------------------------------------------------------------------
  ## Arguments:  p + q == 1

    ##__________ for Rmpfr mpfr() numbers, need a more accurate stirlerr()
    ##--- otherwise cannot get more than double prec.! >>> ./dgamma.R  (using 'DPQmpfr::stirlerrM()')

    stopifnot(is.logical(log), length(log) == 1)
    ## Recycle to common length
    M <- max(length(x), length(n), length(p), length(q))
    r <- double(M) # result
    if(M == 0) return(r)
    ## M >= 1 :
    x <- rep_len(x,  M)
    n <- rep_len(n,  M)
    p <- rep_len(p,  M)
    q <- rep_len(q,  M)

    if(any(B <- p == 0)) r[B] <- ifelse(x[B] ==  0  , .D_1(log), .D_0(log))
    if(any(B <- q == 0)) r[B] <- ifelse(x[B] == n[B], .D_1(log), .D_0(log))
    BB <- !B & (p != 0)
    if(verbose) cat(sprintf("dbinom_raw(): #(bndr, rglr): (%d,%d)\n",
                            sum(!BB), sum(BB)))
    verb1 <- pmax(0L, verbose - 1L)

    if(any(B <- BB & x == 0)) {
        ii <- which(B)
        if(verbose) cat(sprintf("x=0 for i = %s\n", deparse(ii, control="S")))
	if(any(i0 <- n[ii] == 0)) {
            r[ii[i0]] <- .D_1(log)
            ii <- ii[!i0] # for the rest of this clause
        }
        n. <- n[ii]
	lc <- ifelse(p[ii] < 0.1,
                     -bd0(n., n.*q[ii], verbose=verb1) - n.*p[ii],
                     n.*log(q[ii]))

        r[ii] <- .D_exp(lc, log)
        BB <- BB & !B
    }
    if(any(B <- BB & x == n)) {
        ii <- which(B)
        if(verbose) cat(sprintf("x=n for i = %s\n", deparse(ii, control="S")))
        n. <- n[ii]
	lc <- ifelse(q[ii] < 0.1,
                     -bd0(n.,n.*p[ii], verbose=verb1) - n.*q[ii],
                     n.*log(p[ii]))

        r[ii] <- .D_exp(lc, log)
        BB <- BB & !B
    }
    if(any(B <- BB & x < 0 | x > n)) {
        if(verbose) cat(sprintf("have %d x values outside [0,n]\n", sum(B)))
        r[B] <- .D_0(log)
        BB <- BB & !B
    }

    if(any(BB)) {
        ii <- which(BB)
        n <- n[ii]
        x <- x[ii]
        if(verbose) cat(sprintf("main, i = %s\n", deparse(ii, control="S")))
        ##  n*p or n*q can underflow to zero if n and p or q are small.  This
        ##  used to occur in dbeta, and gives NaN as from R 2.3.0.
        lc <- { stirlerr(n, verbose=verb1) - stirlerr(x, verbose=verb1) - stirlerr(n-x, verbose=verb1) -
                    bd0( x , n*p[ii], verbose=verb1) - bd0(n-x, n*q[ii], verbose=verb1)
        }

        ## f = (M.2PI*x*(n-x))/n; could overflow or underflow */
        ##lf <- log(2*pi) + log(x) + log(n-x) - log(n)
        ##                          ---------------- = log((n-x)/n)=log(1 - x/n)
        lf  <- log(2*pi) + log(x) + log1p(-x/n)

        if(verbose) cat(sprintf("  lc=%g, lf=%g ==> lc - 0.5*lf = %g\n", lc,lf,lc - 0.5*lf))

        r[BB] <- .D_exp(lc - 0.5*lf, log)
    }
    r
}

dnbinomR <- function (x, size, prob, log = FALSE, eps = 1e-10)
{
  ## Purpose: R version of R'C level dnbinom() in .../R/src/nmath/dnbinom.c

    x <- floor(x + 1e-7)
    stopifnot(is.logical(log), length(log) == 1,
              0 < prob, prob <= 1, size >= 0, x >= 0)
    M <- max(length(x), length(size), length(prob))
    r <- double(M)
    if(M == 0) return(r)
    x    <- rep_len(x,     M)
    size <- rep_len(size,  M)
    prob <- rep_len(prob,  M)

    ## This is  ** in addition ** to the C code .. and is part of x < eps * size below
    if(any(i0 <- x == 0)) {
        if(any(i0s0 <- i0 & size == 0)) {
            r[i0s0] <- .D_1(log)
        }
        if(any(i0P <- i0 & size > 0)) {  ## x = 0, size > 0
            ## pr(x,...) = pr^n :
            r[i0P] <- if(log) size[i0P]*log(prob[i0P]) else prob[i0P] ^ size[i0P]
        }
           x <-    x[!i0]
        size <- size[!i0]
        prob <- prob[!i0]
    }
    if(any(B <- x < eps * size)) { ## don't use dbinom_raw() but MM's formula
        i <- which(B)
        x. <-    x[i]
        n. <- size[i]
        pr <- prob[i]
        r[!i0][i] <- .D_exp(n. * log(pr) + x. * (log(n.) + log1p(-pr))
                            -lgamma1p(x.) + log1p(x.*(x.-1)/(2*n.)),
                            log)
    }
    if(any(!B)) {
        i <- which(!B)
        x. <- x   [i]
        n. <- size[i]
        pr <- prob[i]
        ans <- dbinom_raw(x = n., n = x.+n., p = pr, q = 1-pr, log = log)
        ## p <- n./(n.+x) ## == 1 if |x| << n.;
        ## better in log case: log(n/(n+x)) = log(1 - x/(n+x))
        r[!i0][i] <- if(log) log1p(-x/(n.+x)) + ans  else  n./(n.+x) * ans
    }
    r
}

dnbinom.mu <- function(x, size, mu, log = FALSE, eps = 1e-10)
{
  ## Purpose: R version of dnbinom_mu() { in .../R/src/nmath/dnbinom.c }
  ## ----------------------------------------------------------------------
  ## Arguments: as  dbinom_mu()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 Jul 2008, 21:54

    x <- floor(x + 1e-7)
    stopifnot(is.logical(log), length(log) == 1,
              mu >= 0, size >= 0, x >= 0)
    M <- max(length(x), length(size), length(mu))
    r <- double(M)
    if(M == 0) return(r)
    x    <- rep_len(x,    M)
    size <- rep_len(size, M)
    mu   <- rep_len(mu,   M)

    if(any(i0 <- x == 0)) {
        i <- which(i0)
        x.  <-   x[i]
        mu. <-  mu[i]
        n. <- size[i]
        ## pr(x,...) = p^n = (n / (n+mu))^n  -- carefully evaluated:
	r[i] <- .D_exp(n. * ifelse(n. < mu.,
                                    log(n./(n.+mu.)),
                                    log1p( - mu./(n.+mu.))),
			 log)
	size <- size[!i0]
	   x <-	   x[!i0]
	  mu <-	  mu[!i0]
    }
    if(length(i <- which(B <- x < eps * size))) { ## don't use dbinom_raw() but MM's formula
        ## log p__r =
        x.  <-   x[i]
        mu. <-  mu[i]
        n. <- size[i]
        r[!i0][i] <- .D_exp(x. * log(n.*mu./(n.+mu.)) - mu.
                            -lgamma1p(x.) + log1p(x.*(x.-1)/(2*n.)),
                            log)
    }
    if(length(i <- which(!B))) {
        x    <-    x[i]
        mu   <-   mu[i]
        size <- size[i]
        ans <- dbinom_raw(x= size, n= x+size,
                          p= size/(size+mu), q= mu/(size+mu),
                          log = log)
        ## p <- size/(size+x) ## == 1 if  |x| << size : can be better in log case:
        ## log(n/(n+x)) = log(1 - x/(n+x))
        r[!i0][i] <- if(log) log1p(-x/(size+x)) + ans else  size/(size+x) * ans
    }
    r
}
