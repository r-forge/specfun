#### Vectorized DCDFLIB distributions
####
#### t, F, Chi^2  --- central & noncentral


### Original was :
## 12/12/97
## Splus functions to use DCDFLIB cdftnc
## Donald H. MacQueen
## macqueen1@llnl.gov

is1 <- function(x) abs(x-1) <= 100*.Machine$double.eps

ptnc <- function(t,df=stop("no df arg"),ncp=0, lower.tail = TRUE)
    tncV(which=1,t=t,df=df,ncp=ncp)[[c("q","p")[1+as.logical(lower.tail)]]]

qtnc <- function(p,df=stop("no df arg"),ncp=0, lower.tail = TRUE) {
    if(lower.tail) tncV(which=2,p=p,   q=1-p, df=df,ncp=ncp)$t
    else           tncV(which=2,p=1-q, q=q,   df=df,ncp=ncp)$t
}


## This is further generalized in 'tfchiV'  --> BELOW
tncV <- function(which = 1, p= 1-q, q= 1-p,
                 t,df, ncp=0,
                 status=0,bound=0)
{

    which <- as.integer(which[1])
    if(which < 1 || which > 4) stop("'which' must be in  {1,2,3,4}")
    switch(which,
       {p <- q <- 0}, ## 1 : calculate p and q
       { t <- 0    }, ## 2 : calculate t (the inverse cdf)
       { df <- 0   }, ## 3 : calculate df
       { ncp <- 0  }  ## 4 : calculate ncp
           )
    if(which!=1 && (( missing(p) && missing(q)) ||
                    (!missing(p) && !missing(q) && any(!is1(p+q)))))
        stop("must specify either 'p' or 'q'; if both, p+q = 1")
    if(which!=3 && (missing(df) || !is.numeric(df) || any(df <= 0)))
        stop("'df' must be > 0")
    if(which!=2 && (missing(t) || !is.numeric(t)))
        stop("'t' must be specified (numeric).")

    lens <- c( length(p),length(t),length(df),length(ncp) )
    len <- max(lens)
    if (any(is.na(match(lens, c(1,len)))))
        warning("tncV(): lengths of arguments not all equal or 1\n")

    vt <- .C(C_V_cdftnc,
             which = which,
             p = rep(as.double(p),  length=len),
             q = rep(as.double(q),  length=len),
             t = rep(as.double(t),  length=len),
             df= rep(as.double(df), length=len),
             pnonc = rep(as.double(ncp),length=len),
             status= rep(as.integer(status),length=len),
             bound = rep(as.double(bound),length=len),
             len = as.integer(len))
    if (any(vt$status != 0)) {
        cat("Warning: error flag returned by DCDFLIB cdftnc\n")
        cat("Call function tncV() with which =", which,
            " and same values for t, df, and ncp\n")
        cat("status\n", vt$status, "\n")
        cat("bound\n", vt$bound, "\n")
        vt <- NULL
    }
    vt
} ## end tncV

tfchiV <- function(which=1,
                   dist = c("t",   "F",  "chi2",
                            "tnc", "Fnc","chi2nc"),
                   p= 1-q, q= 1-p,
                   t,df, df2=0, ncp=0,
                   status=0, bound=0)
{
    ## Purpose: Vectorized (noncentral) t-, F-, or chi^2-  distribution, etc.
    ##  Goal / Approach / Advantages:
    ##	1) manual (and transparent!) vectorization
    ##    2) As per design-goal of DCDFLIB:  Specify all but one "parameter"
    ##         --> "solve" for the missing one !
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 10 Apr 99

    dists <- as.character(as.list(formals()[["dist"]])[-1])# defaults above
    if(length(dist) != 1 || !is.character(dist))
        stop(paste("'dist' must be character, ",
                   "(possibly an abbreviation of) one of\n \"",
                   paste(dists, collapse='", "'),'"',sep=""))
    dist <- match.arg(dist)

    ## do we have a NON-central at all?
    nc <- nchar(dist)
    non.cent <- substr(dist,nc-1, nc) == "nc"

    ## UNFINISHED ....

    which <- as.integer(which[1])
    if(which < 1 || which > 4) stop("'which' must be in  {1,2,3,4}")
    if(which == 4 && !non.cent)
        stop(paste("which=4  does not make sense for", dist))
    switch(which,
       {p <- q <- 0}, ## 1 : calculate p and q
       { t <- 0    }, ## 2 : calculate t (the inverse cdf)
       { df <- 0   }, ## 3 : calculate df  ( & df2 ????)
       { ncp <- 0  }  ## 4 : calculate ncp
           )
    if(which!=1 && (( missing(p) && missing(q)) ||
                    (!missing(p) && !missing(q) && any(!is1(p+q)))))
        stop("must specify either 'p' or 'q'; if both, p+q = 1")
    if(which!=3 && (missing(df) || !is.numeric(df) || any(df <= 0)))
        stop("'df' must be > 0")
    if(which!=2 && (missing(t) || !is.numeric(t)))
        stop("'t' must be specified (numeric).")

    lens <- c( length(p),length(t),length(df),length(df2),length(ncp) )
    len <- max(lens)
    if(any(is.na(match(lens, c(1,len)))))
        warning("tfchiV(*, method=", dist,"): lengths of arguments not all equal or 1\n")

    p <- rep(as.double(p),  length=len)
    q <- rep(as.double(q),  length=len)
    t <- rep(as.double(t),  length=len)
    df<- rep(as.double(df), length=len)
    if(non.cent) pnonc <- rep(as.double(ncp),length=len)
    status<- rep(as.integer(status),length=len)
    bound <- rep(as.double(bound),length=len)
    len <- as.integer(len)
    vt <- switch(dist,
           t   = { .C(C_V_cdft, which, p=p, q=q, t=t, df=df,
                      status=status, bound = bound, len) },
           F   = { .C(C_V_cdff, which,p=p,q=q,t=t,df1=df,df2=df2,
                      status=status, bound = bound, len) },
           chi2= { .C(C_V_cdfchi, which, p=p, q=q, t=t, df=df,
                      status=status, bound = bound, len) },
           tnc = { .C(C_V_cdftnc, which, p=p, q=q, t=t, df=df,
                      pnonc=pnonc, status=status, bound = bound, len) },
           Fnc = { .C(C_V_cdffnc, which, p=p, q=q, t=t, df1=df,df2=df2,
                      pnonc=pnonc, status=status, bound = bound, len) },
           chi2nc={ .C(C_V_cdfchn, which, p=p, q=q, t=t, df=df,
                       pnonc=pnonc, status=status, bound = bound, len) }
           )

    if (any(vt$status != 0)) {
        warning("error flag returned by DCDFLIB's cdf<", dist,">()\n", sep="")
        cat("called with which =", which,
            " and .?.same.?. values for t, df, and ncp\n")
        cat("status\n", vt$status, "\n")
        cat("bound\n", vt$bound, "\n")
        #vt <- NULL
    }
    vt
}
