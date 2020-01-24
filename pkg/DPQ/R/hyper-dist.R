#### Hypergeometric Distribution
#### ===========================
#### R ( / S) functions :
####	dhyper(x, m, n, k)
####	phyper(q, m, n, k)
###
### Parametrization / Nomenclature :
###
### JohKK  R/S        C    AS152
### ~~~~~  ~~~       ~~
###  Np   = m     =  NR   = nn   : # { WHITE (Red) balls   in the urn }
###  N-Np = n     =  NB   = mm-nn: # { BLACK       balls   in the urn }
###   n   = k     =   n   = kk   : # { balls drawn from the urn }
###         x     =   x   = ll   : # { WHITE (Red) balls AMONG the `n' drawn}
###  ----  ---      -----
###   N   =m+n    = NR+NB = mm   : the TOTAL # { balls   in the urn }
###   p   =m/(m+n)= NR/N  = nn/mm: probability of WHITE(Red)

### where
###  JohKK := Johnson, Kotz & Kemp (1992)
###           Univariate Discrete Distributions, 2nd Ed.
###           [Chapter 6 -- Hypergeometric Distributions]
###  R/S   := the [dpqr]hyper() functions defined in R {and already S}

### Some tests (no longer problematic):
###  ./dpq-functions/Knuesel-splus-ex.R  [ -->  phyper() ]

#### Used to be part of /u/maechler/R/MM/NUMERICS/hyper-dist.R -- till Jan.24 2020

#### First version committed is the code "as in Apr 1999":
####  file | 11782 | Apr 22 1999 | hyper-dist.R

phyper.appr.as152 <- function(q, m, n, k)
{
  ## Purpose: Normal Approximation to cumulative Hyperbolic
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 19 Apr 99, 23:39
  kk <- n
  nn <- m
  mm <- m+n
  ll <- q
  mean <- kk * nn / mm
  sig <- sqrt(mean * (mm - nn) / mm * (mm - kk)/(mm - 1) )
  pnorm(ll + 1/2, mean=mean, sd = sig)
}

phyper.Ibeta <- function(q, m, n, k)
{
  ## Purpose: Pearson [Incompl.Beta] Approximation to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.90), p.260 -- Bol'shev (1964)
  ## ----------------------------------------------------------------------

  Np <- m; N <- n + m; n <- k; x <- q
  p <- Np/N ; np <- n*p
  xi <- (n + Np -1 - 2*np) / (N-2)

  d.c <- (N-n)*(1-p) + np - 1
  cc <- n*(n-1)*p*(Np-1) / ((N-1)*d.c)
  lam <- (N-2)^2 *np *(N-n)*(1-p) / ((N-1)*d.c*(n + Np-1 - 2*np))
  pbeta(1 - xi,
        lam - x + cc,
        x - cc + 1)
}

phyper.1molenaar <- function(q, m, n, k)
{
  ## Purpose: Normal Approximation to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.91), p.261
  ## ----------------------------------------------------------------------

  Np <- m; N <- n + m; n <- k; x <- q
  pnorm(2/sqrt(N-1) *
        (sqrt((x + 1)*(N - Np - n + x + 1)) -
         sqrt((n - x)*(Np - x)) ))
}

phyper.2molenaar <- function(q, m, n, k)
{
  ## Purpose: Normal Approximation to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.92), p.261
  ## ----------------------------------------------------------------------
  Np <- m; N <- n + m; n <- k; x <- q
  pnorm(2/sqrt(N) *
        (sqrt((x + .75)*(N - Np - n + x + .75)) -
         sqrt((n - x -.25)*(Np - x - .25)) ))
}
phyper.Peizer <- function(q, m, n, k)
{
  ## Purpose: Peizer's extremely good Normal Approx. to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.93) & (6.94), p.261 __CORRECTED__
  ## ----------------------------------------------------------------------
  Np <- m; N <- n + m; n <- k; x <- q
  ## (6.94) -- in proper order!
  nn <- Np			;  n. <- Np     + 1/6
  mm <- N - Np                  ;  m. <- N - Np + 1/6
  r <- n                        ;  r. <- n      + 1/6
  s <- N - n                    ;  s. <- N - n  + 1/6
                                   N. <- N      - 1/6
  A <- x + 1/2                  ;  A. <- x      + 2/3
  B <- Np - x - 1/2             ;  B. <- Np - x - 1/3
  C <- n  - x - 1/2             ;  C. <- n  - x - 1/3
  D <- N - Np - n + x + 1/2     ;  D. <- N - Np - n + x + 2/3

  n <- nn
  m <- mm
  ## After (6.93):
  L <-
    A * log((A*N)/(n*r)) +
    B * log((B*N)/(n*s)) +
    C * log((C*N)/(m*r)) +
    D * log((D*N)/(m*s))
  ## (6.93) :
  pnorm((A.*D. - B.*C.) / abs(A*D - B*C) *
        sqrt(2*L* (m* n* r* s* N.)/
                  (m.*n.*r.*s.*N )))
  # The book wrongly has an extra "2*" before `m* ' (after "2*L* (" ) above
}

### Binomial Approximation(s) ==========
### ======== ~~~~~~~~~~~~~~~~ ==========

##  JohKK  R     C
##  ~~~~~  ~~    ~~
##  Np   = m  =  NR  : # { WHITE (Red) balls   in the urn }
##  N-Np = n  =  NB  : # { BLACK       balls   in the urn }
##   n   = k  =   n  : # { balls drawn from the urn }
##         x  =   x  : # { WHITE (Red) balls AMONG the `n' drawn}

##' Molenaar(1970a) -- formula (6.80) in JohKK -- improved p for Bin(n,p) approximation:
hyper2binom.p <- function(x, m,n,k) {
    N <- m+n
    p <- m / N # "Np / N"
    N.n <- N - (k-1)/2
    ## p^{+} =
    (m - x/2)/N.n - k*(x - k*p - 1/2) / (6 * N.n^2)
}

### NOTA BENE/TODO: (this paragraph has not __yet__ been applied below at all !!)
### ---------
### JohKK, p.258 (top) mention the *four* different binomial approximations
### for a given hypergeometric, and then
### "Brunk et al. (1968) .. support the opinion of ...(1961) that it is
### best to use the binomial with smallest power parameter, that is,
### n' = min(n, Np, N-Np, N-n)
### --------------------------
###   ( x ; (  n , p= Np/N) )  <--> ( x  ; (Np , p= n /N)) <-->
###   (n-x; (N-Np, p= n /N) )  <--> (Np-x; (N-n, p= Np/N))
### translated to R's notation:
###   ( q ; ( k ,  p= m/(m+n)) <--> ( q  ; (m,     p= k/(m+n))) <-->
###   (k-q; ( n ,  p= k/(m+n)) <--> (m-q ; (m+n-k, p= m/(m+n)))

##' The support of the hypergeometric distrib. as a function of its parameters
q.mnk.hyper <- function(m,n,k) max(0, k-n) : min(k, m)

##' The two symmetries <---> the four different ways to compute :
phypers <- function(m,n,k, q = q.mnk.hyper(m,n,k)) {
    N <- m+n
    pm <- cbind(ph = phyper(q,     m,  n , k), # 1 = orig.
                p2 = phyper(q,     k, N-k, m), # swap m <-> k (keep N = m+n)
                ## "lower.tail = FALSE"  <==>  1 - p..(..)
                Ip2= phyper(m-1-q, N-k, k, m, lower.tail=FALSE),
                Ip1= phyper(k-1-q, n,   m, k, lower.tail=FALSE))

    ## check that all are (approximately) the same :
    stopifnot(all.equal(pm[,1], pm[,2]),
              all.equal(pm[,2], pm[,3]),
              all.equal(pm[,3], pm[,4]))
    list(q = q, phyp = pm)
}

phyper.bin.Molenaar <-
phyper.bin.Molenaar.1 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(q, size = k, prob = hyper2binom.p(q, m,n,k),
           lower.tail=lower.tail, log.p=log.p)
phyper.bin.Molenaar.2 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    ## swap k ('n') with m ('Np') -- but with R's notation n=N-Np changes too:
    pbinom(q, size = m, prob = hyper2binom.p(q, k, n-k+m, m),
           lower.tail=lower.tail, log.p=log.p)
phyper.bin.Molenaar.3 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE) {
    ## "Ip2"
    pbinom(m-1-q, size = m, prob = hyper2binom.p(m-1-q, m+n-k, k, m),
           lower.tail = !lower.tail, log.p=log.p)
    ##                 ===
}
phyper.bin.Molenaar.4 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE) {
    ## "Ip1"
    pbinom(k-1-q, size = k, prob = hyper2binom.p(k-1-q, n, m, k),
           lower.tail = !lower.tail, log.p=log.p)
    ##                 ===
}

## Now, for completeness, also the simple binomial approximations:
phyper.bin.1 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(q, size = k, prob = m/(m+n), lower.tail=lower.tail, log.p=log.p)
phyper.bin.2 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    ## swap k ('n') with m ('Np') -- but with R's notation n=N-Np changes too:
    pbinom(q, size = m, prob = k/(m+n), lower.tail=lower.tail, log.p=log.p)
phyper.bin.3 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(m-1-q, size = m, prob = (m+n-k)/(m+n), lower.tail = !lower.tail, log.p=log.p)
phyper.bin.4 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(k-1-q, size = k, prob = n/(m+n), lower.tail = !lower.tail, log.p=log.p)


phyper.bin.allM <- function(m, n, k, q = q.mnk.hyper(m,n,k),
                            lower.tail=TRUE, log.p=FALSE)
{
    cbind(pM1 = phyper.bin.Molenaar.1(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM2 = phyper.bin.Molenaar.2(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM3 = phyper.bin.Molenaar.3(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM4 = phyper.bin.Molenaar.4(q, m, n, k, lower.tail=lower.tail, log.p=log.p))
}

phyper.bin.all <- function(m, n, k, q = q.mnk.hyper(m,n,k),
                           lower.tail=TRUE, log.p=FALSE)
{
    cbind(
          p1 = phyper.bin.1(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          p2 = phyper.bin.2(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          p3 = phyper.bin.3(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          p4 = phyper.bin.4(q, m, n, k, lower.tail=lower.tail, log.p=log.p),

          pM1 = phyper.bin.Molenaar.1(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM2 = phyper.bin.Molenaar.2(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM3 = phyper.bin.Molenaar.3(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM4 = phyper.bin.Molenaar.4(q, m, n, k, lower.tail=lower.tail, log.p=log.p))

}

dhyper.bin.Molenaar <- function(x, m, n, k, log=FALSE)
    dbinom(x, size=k, prob = hyper2binom.p(x, m,n,k), log=log)

### ----------- ----------- lfastchoose() etc -----------------------------

lfastchoose <- function(n,k) lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

my.lchoose <- function(n,k) {
  h <- 0.5; n <- floor(n+h); k <- floor(k+h); lfastchoose(n,k)
}

## in  math/gamma.c :
p1 <- c(0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2)
## (nowhere used ??)


## The  Bernoulli Numbers:
Bern <- function(n, verbose = getOption("verbose", FALSE))
{
  ## Purpose: n-th Bernoulli number -- exercise in cashing
  ## ----------------------------------------------------------------------
  ## Arguments: n >= 0   B0 = 1, B1 = -1/2,  B2 = 1/6,  B4 = -1/30,...
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 26 Apr 97, 18:17
  n <- as.integer(n)
  if(n < 0) stop("'n'  must not be negative")
  if(n== 0) 1 else
  if(n==1) -1/2 else
  if(n %% 2 == 1) 0 else { ##-- maybe use 'cache': .Bernoulli[n] == Bern(2 * n)
    n2 <- n %/% 2
    do.new <- !exists(".Bernoulli", mode='numeric')
    if(do.new) .Bernoulli <<- numeric(0)
    if(do.new || length(.Bernoulli) < n2) { ##-- Compute  Bernoulli(n)
      if(verbose) cat("n=",n,": computing", sep='', "\n")
      Bk <- k0 <- seq(length=n2-1)
      if(n2 > 1) {
	for(k in k0) Bk[k] <- Bern(2*k)
	k0 <- 2 * k0
      }
      .Bernoulli[n2] <<- - sum( choose(n+1, c(0,1,k0)) * c(1,-1/2, Bk)) / (n+1)
    }
    .Bernoulli[n2]
  }
}

### as.lgamma() --- Asymptotic log gamma function :

as.lgamma <- function(x, n)
{
  ## Purpose: asymptotic log gamma function
  ## ----------------------------------------------------------------------
  ## Arguments: x >~ 3 ; n: number of terms in sum
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 27 Apr 97, 22:18
  s <- (x-1/2)*log(x) -x + log(2*pi)/2
  if(n>=1) {
    Ix2 <- 1/(x*x)
    Bern(2*n) #-> assigning  .Bernoulli[1:n]
    k <- 1:n
    Bf <- rev(.Bernoulli[k] / (2*k*(2*k-1)))
    bsum <- Bf[1]
    for(i in k[-1])
      bsum <- Bf[i] + bsum*Ix2
    s + bsum/x
  } else s
}


### R version of  {old} C version in  <R>/src/nmath/phyper.c :

## /usr/local/app/R/R-MM-release/src/nmath/phyper.c

phyperR <- function(x, NR, NB, n)
{
    x <- floor(x)
    NR <- floor(NR + 0.5)
    NB <- floor(NB + 0.5)
    N <- NR + NB
    n <- floor(n + 0.5)
    if (NR < 0 || NB < 0 || n < 0 || n > N) {
        stop("domain ERROR in phyper")
    }
    xstart <- max(0, n - NB)
    xend <- min(n, NR)
    if(x < xstart) return(0.0)
    if(x >= xend) return(1.0)
    xr <- xstart
    xb <- n - xr
    ltrm <- lchoose(NR, xr) + lchoose(NB, xb) - lchoose(N, n)
    ##o term <- exp(ltrm)
    NR <- NR - xr
    NB <- NB - xb
    Sm <- s2 <- 0.0
    while(xr <= x) {
        ##o Sm <- Sm + term
        s2  <- s2  + exp(ltrm)
        xr <- xr+1
        NB <- NB+1
        ff <- NR * xb / (xr * NB) ## (NR / xr) * (xb / NB)
        ltrm <- ltrm + log(ff)
        ##o term <- term * ff
        xb <- xb-1
        NR <- NR-1
    }
    s2 ##o return(Sm,s2)
}
