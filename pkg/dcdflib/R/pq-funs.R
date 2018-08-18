####--- Provide all the  p<dist> and q<dist> functions

## cdfbet(WHICH,p,q, x,y, a,	b,	STATUS,bound) c.d.f. Beta
## cdfbin(WHICH,p,q, s,xn, pr,	ompr,	STATUS,bound) c.d.f. Binomial
## cdfchi(WHICH,p,q, x,df,		STATUS,bound) c.d.f. Chisquare
## cdfchn(WHICH,p,q, x,df, pnonc,	STATUS,bound) c.d.f. Non-c.Chisquare
## cdff  (WHICH,p,q, f,dfn, dfd,	STATUS,bound) c.d.f. F
## cdffnc(WHICH,p,q, f,dfn, dfd, pnonc, STATUS,bound) c.d.f. Non-central F
## cdfgam(WHICH,p,q, x,shape, scale,	STATUS,bound) c.d.f. Gamma
## cdfnbn(WHICH,p,q, s,xn, pr,	ompr,	STATUS,bound) c.d.f. Negative Binomial
## cdfnor(WHICH,p,q, x,mean, sd,	STATUS,bound) c.d.f. Normal
## cdfpoi(WHICH,p,q, s,xlam, 		STATUS,bound) c.d.f. Poisson
## cdft  (WHICH,p,q, t,df, 		STATUS,bound) c.d.f. T
## cdftnc(WHICH,p,q, t,df, pnonc,	STATUS,bound) c.d.f. Non-central T

##-- New: use cdfbetV()  etc for  VECTORIZED operation !

##--- There's no NON-centrality available here
##--- Use  non-central F  for this !
pbetaB <- function(x = 1-y, a, b, y = 1-x, lower.tail = TRUE)
{
  if(any(a<=0) || any(b<=0)) stop("'a' and 'b' must both be positive")
  if(any(!is1(x+y))) stop("'x == 1-y' is required")

  res <- cdfbetV(which=1, x=x,y=y,a=a,b=b,
                 p=0,q=1,status=0,bound=0)[c("p","q","status","bound")]
  if(any(Prob <- res$status != 0)) {
      cc <- function(...) paste("(", paste(..., collapse =", "),")", sep="")
      warning(paste("pbetaB(.): status  for x[",cc(which(Prob)),"] =",
                    cc(res$status[Prob]),
                    "; bound=", cc(res$bound[Prob]), sep=""))
  }
  res[[c("q","p")[1+as.logical(lower.tail)]]]
}

qbetaB <- function(p, a, b, lower.tail = TRUE) {
  if(any(a<=0) || any(b<=0)) stop("'a' and 'b' must both be positive")
  if(lower.tail) q <- 1-p else {
    q <- p
    p <- 1-q
  }
  res <- cdfbetV(which = 2, x=0,y=1, a=a,b=b, p=p, q=q,
                 status=0,bound=0)[c("x","y","status","bound")]
  if(any(Prob <- res$status != 0)) {
      cc <- function(...) paste("(", paste(..., collapse =", "),")", sep="")
      warning(paste("qbetaB(.): status for x[",cc(which(Prob)),"] =",
                    cc(res$status[Prob]),
                  "; bound=", cc(res$bound[Prob]), sep=""))
  }
  res$x # res$y  might be interesting as well ..
}

pgammaB <- function(x, shape, scale=1, lower.tail = TRUE)
{
    if(any(shape<=0) || any(scale<=0))
        stop("'shape' and 'scale' must both be positive")

    res <- cdfgamV(which=1, p=0,q=1, x=x, shape=shape,scale=scale,
                   status=0,bound=0)[c("p","q","status","bound")]
    if(any(Prob <- res$status != 0)) {
        cc <- function(...) paste("(", paste(..., collapse =", "),")", sep="")
        warning(paste("pgammaB(.): status for x[",cc(which(Prob)),"] =",
                      cc(res$status[Prob]),
                      "; bound=", cc(res$bound[Prob]), sep=""))
    }
    res[[c("q","p")[1+as.logical(lower.tail)]]]
}

qgammaB <- function(p, shape, scale=1, q=1-p, lower.tail = TRUE) {
    if(any(shape<=0) || any(scale<=0))
        stop("'shape' and 'scale' must both be positive")
    if(lower.tail) q <- 1-p else {
        q <- p
        p <- 1-q
    }
    res <- cdfgamV(which = 2,  p=p, q=q, x=0, shape=shape, scale=scale,
                   status=0,bound=0)[c("x","y","status","bound")]
    if(any(Prob <- res$status != 0)) {
        cc <- function(...) paste("(", paste(..., collapse =", "),")", sep="")
        warning(paste("qgammaB(.) : status for x[",cc(which(Prob)),"] =",
                      cc(res$status[Prob]),
                      "; bound=", cc(res$bound[Prob]), sep=""))
    }
    res$x
}


pbinomB <- function(x, n, p, lower.tail = TRUE) {
}

qbinomB <- function(prob, n, p, q=1-prob, lower.tail = TRUE) {
}


pfB <- function(x, n1, n2, ncp=0, lower.tail = TRUE) {
}
qfB <- function(p, n1, n2, ncp=0, q=1-p, lower.tail = TRUE) {
}



pchisqB <- function(x, df, ncp=0, lower.tail = TRUE) {
}
qchisqB <- function(p, df, ncp=0, q=1-p, lower.tail = TRUE) {
}




pnormB <- function(x, mean=0, sd=1, lower.tail = TRUE) {
}

qnormB <- function(p, mean=0, sd=1, q=1-p, lower.tail = TRUE) {
}


##--> ./nct.R  for the  (non)central t functions
ptB <- function(x, df, ncp=0, lower.tail = TRUE) {
}
qtB <- function(p, df, ncp=0, q=1-p, lower.tail = TRUE)
{}

