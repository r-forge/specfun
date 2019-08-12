#### Approximations of   qnchisq() -*- delete-old-versions: never -*-
####
#### R simulation of C's qnchisq()  <<< not yet >>
####
####  =========================
####  Function definitions only
####  =========================

###-- 1 -- Normal Approximations ---------------------------------------

##- Date: Wed, 19 May 1999 16:47:23 -0400 (EDT)
##- From: Ranjan Maitra <maitra@math.umbc.edu>
##- To: S-news <s-news@wubios.wustl.edu>
##- Subject: Re: [S] quantiles of non-central chisq

##- Many thanks to Jim Stapleton for answering my request for a function to
##- calculate the quantiles of a non-central chi-squared distribution.

## MM: Where is this approx. from?

##- The function is as follows:
qchisqAppr.0 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    df + ncp + z.p * (2 * df + 4 * ncp)^0.5
}

qchisqAppr.1 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## (1+b) * q^chisq(x,df);
    ## better (?) to use q^chisq(x, df*); -> qchisqAppr.2
    (1 + ncp/(ncp +df))* df * (1 - 2/(9*df) + z.p*sqrt(2/(9*df)))^3
}

qchisqAppr.2 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## use = q^chisq(x, df*); df* = a/(1+b); a=(df+ncp); b=ncp/(df+ncp)
    a <- df + ncp
    b1 <- 1+ ncp/a # = 1+b
    df. <- a/b1
    alp <- 2/(9*df.)
    a * (1 - alp + z.p *sqrt(alp))^3
}

qchisqAppr.3 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## use = q^chisq(x, df*); df* = a/(1+b); a=(df+ncp); b=ncp/(df+ncp)
    a <- df + ncp
    df. <- a/(1+ ncp/a)
    alp <- 2/(9*df.)
    df. * (1 - alp + z.p *sqrt(alp))^3
}


## -- Cornish-Fisher (inverse Edgeworth) Expansions:
qchisqApprCF1 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## 1-term Cornish-Fisher Expansion
    mu <- df + ncp
    sig2 <- 2*(df + 2*ncp)
    gam1 <- 8*(df + 3*ncp) * sig2^(-3/2)
    mu + sqrt(sig2)* (z.p + gam1/6 * (z.p^2 - 1))
}

qchisqApprCF2 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## 2-term Cornish-Fisher Expansion
    mu <- df + ncp
    sig2 <- 2*(df + 2*ncp)
    gam1 <- 8*(df + 3*ncp) * sig2^(-3/2)
    gam2 <- 48*(df + 4*ncp) / sig2^2
    mu + sqrt(sig2)* (z.p + gam1/6 * (z.p^2 - 1) +
                      gam2/24* z.p*(z.p^2-3) - gam1^2/36* z.p*(2*z.p^2 -5))
}



###-- 2 -- Chisq (Central) Approximations ---------------------------

qnchisqPatnaik <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Basically Abramowitz & Stegun "nu > 30" formula  (26.4.30)
    ## Also == Patnaik(1949)'s as reported in Johnson,Kotz,Bala.. p.462
    a <- df + ncp
    b1 <- 1 + ncp / a
    nu. <-  a/b1 # nu* = (nu + lam)^2 / (nu + 2*lam)
    b1 * qchisq(p, df = nu., lower.tail= lower.tail, log.p= log.p)
}

qnchisqPearson <- function(p, df, ncp = 0,
                       lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Pearson(1959) approximation to pnchisq() - using pchisq()
    ##  Very good for df > 2--5  and much better at *right tail*
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 27 Feb 2004, 15:37

    n2 <- df + 2*ncp
    n3 <- n2 + ncp ## = df + 3*ncp
    -(ncp*ncp /n3) + (n3/n2)*
        qchisq(p, df = n2^3/(n3*n3), lower.tail=lower.tail, log.p=log.p)
}

qchisqCappr.2 <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Johnson,Kotz,Bala.. p.466  (29.66)
    ## is  *VERY* bad for small df or small df/ncp  (df = 0.1; df/ncp < 1};
    ## but quite strong for small ncp, or small ncp/df , even (1 , 1)
    c2 <- qchisq(p, df, lower.tail= lower.tail, log.p= log.p)
    l.n <- ncp / df # lambda / nu
    c2 * (1 + l.n *(1 + 0.5*l.n*(1 - c2 / (df+2))))
}

## Nice inversions of Normal approximations in ./pnchisq.R :
## Abdel-Aty (1954) .. Biometrika
qnchisqAbdelAty <-  function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Abdel-Aty(1954) "first approx." (Wilson-Hilferty) approximation to pnchisq()
    ## Johnson,Kotz,...: Has it *WRONGLY*  in eq. (29.61a), p.463
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 2018-08-16
    ## Explicit inversion of pnchisqAbdelAty()
    r <- df + ncp # = ν + λ  ( = k + λ  on Wikip.)
    n2 <- r + ncp # = ν + 2λ
    V <- 2*n2 / (9*r^2)
    qNp <- qnorm(p, mean = 1-V, sd = sqrt(V), lower.tail=lower.tail, log.p=log.p)
    r * qNp^3
}

qnchisqSankaran_d <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Sankaran(1959,1963) according to Johnson et al. p.463, (29.61d):
    ## Explicit inversion of pnchisqSankaran_d()
    r <- df + ncp  #              = ν + λ
    n2 <- r + ncp  # = df + 2*ncp = ν + 2λ
    n3 <- n2 + ncp # = df + 3*ncp = ν + 3λ
    q1 <- n2/r^2   # = (n + 2λ)/(n+λ)²
    h1 <- - 2/3 * r*n3/n2^2 # =  h - 1
    h <- 1 + h1
    mu <- 1 + h*h1*q1*(1 + (h-2)*(1-3*h)*q1/2)
    V  <- 2*h^2 * q1 *(1 +  h1 * (1-3*h)*q1)
    qNp <- qnorm(p, mean = mu, sd = sqrt(V), lower.tail=lower.tail, log.p=log.p)
    r * qNp^(1/h)
}



###-- 3 -- Inversion Algorithms ---------------------------------------


##-- The next ones all build on using  Newton() on pchisq() / dchisq();

## This function finds the quantile corresponding to p
## for the noncentral-chisquare distribution with
##  df degrees of freedom and noncentrality parameter delta = ncp.
qchisq2 <- function(p, df, ncp = 0, eps = 0.0001)
{
    newton(x0 = qchisqAppr.0(p, df, ncp),
           G = function(x, z) pchisq(x, df = z[2], ncp = z[3]) - z[1],
           g = function(x, z) dchisq(x, df = z[2], ncp = z[3]),
           eps = eps, z = c(p, df, ncp))[1]
}

qchisqA <- function(p, df, ncp = 0, eps = 0.0001, qIni = qchisqAppr.0)
{
    ## same as qchisq2() but returning  "all"
    newton(x0 = qIni(p, df=df, ncp=ncp),
           G = function(x, z) pchisq(x, df = z[2], ncp = z[3]) - z[1],
           g = function(x, z) dchisq(x, df = z[2], ncp = z[3]),
           z = c(p, df, ncp), eps = eps, max.iter = 100, give.all = TRUE)
}

newton <- function(x0, G, g, z = 0, eps = 0.001, max.dx = 1, max.iter = 20,
                   give.all = FALSE)
{
    ## Given the function G and its derivative g,
    ## newton uses the Newton method, beginning at x0,
    ## to find a point xp at which G is zero. G and g
    ## may each depend on a parameter z. The result is the
    ##  3-component vector of an approximation x0 of xp,
    ##  G(x0, z), and the number ("iter") of max.iter.
    ##  x0 satisfies  abs(G(x0, z)) < eps.

    ##  give.all = TRUE  to also get  the vectors of consecutive values of
    ##  x0 and G(x0, z);

    iter <- 0
    Gx <- G(x0, z)
    if(give.all) {
        x0vec <- x0
        Gxvec <- Gx
    }
    conv <- TRUE
##-     while(abs(Gx) > eps) {
##-         r <- Gx/g(x0, z)
    while(abs(r <- Gx/g(x0, z)) > eps) {
        d <- min(abs(r), max.dx) * sign(r)
        x0 <- x0 - d
        Gx <- G(x0, z)
        if(give.all) {
            Gxvec <- c(Gxvec, Gx)
            x0vec <- c(x0vec, x0)
        }
        if((iter <- iter + 1) > max.iter) {
            warning("iter >", format(max.iter)) ; conv <- FALSE; break
        }
    }
    if(give.all)
        list(x = x0, x.vec = x0vec, G.vec = Gxvec, it = iter, converged = conv)
    else c(x = x0, G = Gx, it =iter, converged = conv)
} ## newton()
