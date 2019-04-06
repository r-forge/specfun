### From copula package [not yet exported there]. as of Dec.2014:

Re_ <- function(z) if(is.complex(z)) Re(z) else z

##' Properly compute log(x_1 + .. + x_n) for a given (n x d)-matrix of n row
##' vectors log(x_1),..,log(x_n) (each of dimension d)
##' Here, x_i > 0  for all i
##' @title Properly compute the logarithm of a sum
##' @param lx (n,d)-matrix containing the row vectors log(x_1),..,log(x_n)
##'        each of dimension d --- can be "mpfrArray" !!!!
##' @param l.off the offset to substract and re-add; ideally in the order of
##'        the maximum of each column
##' @return log(x_1 + .. + x_n) [i.e., OF DIMENSION d!!!] computed via
##'         log(sum(x)) = log(sum(exp(log(x))))
##'         = log(exp(log(x_max))*sum(exp(log(x)-log(x_max))))
##'         = log(x_max) + log(sum(exp(log(x)-log(x_max)))))
##'         = lx.max + log(sum(exp(lx-lx.max)))
##'         => VECTOR OF DIMENSION d
##' @author Marius Hofert, Martin Maechler
lsum <- function(lx, l.off = apply(Re_(lx), 2, max)) {
    ## do not use cbind or rbind here, since it is not clear if the user specified
    ## only one vector log(x) or several vectors of dimension 1 !!!
    stopifnot(length(dim(lx)) == 2L) # is.matrix(.) generalized
    l.off + log(colSums(exp(lx - rep(l.off, each=nrow(lx)))))
}

##' Properly compute log(x_1 + .. + x_n) for a given matrix of column vectors
##' log(|x_1|),.., log(|x_n|) and corresponding signs sign(x_1),.., sign(x_n)
##' Here, x_i is of arbitrary sign
##' @title compute logarithm of a sum with signed large coefficients
##' @param lxabs (d,n)-matrix containing the column vectors log(|x_1|),..,log(|x_n|)
##'        each of dimension d
##' @param signs corresponding matrix of signs sign(x_1), .., sign(x_n)
##' @param l.off the offset to substract and re-add; ideally in the order of max(.)
##' @param strict logical indicating if it should stop on some negative sums
##' @return log(x_1 + .. + x_n) [i.e., of dimension d] computed via
##'         log(sum(x)) = log(sum(sign(x)*|x|)) = log(sum(sign(x)*exp(log(|x|))))
##'         = log(exp(log(x0))*sum(signs*exp(log(|x|)-log(x0))))
##'         = log(x0) + log(sum(signs* exp(log(|x|)-log(x0))))
##'         = l.off   + log(sum(signs* exp(lxabs -  l.off  )))
##' @author Marius Hofert and Martin Maechler
lssum <- function (lxabs, signs, l.off = apply(Re_(lxabs), 2, max), strict = TRUE) {
    stopifnot(length(dim(lxabs)) == 2L) # is.matrix(.) generalized
    sum. <- colSums(signs * exp(lxabs - rep(l.off, each=nrow(lxabs))))
    if(anyNA(sum.) || any(sum. <= 0))
        (if(strict) stop else warning)("lssum found non-positive sums")
    l.off + log(sum.)
}


##' Compute sqrt(1 + z^2) in a numerical stable way,
##' notably for |z| << 1  (|z| := Arg(z) = ph z  for complex z)
##'
##' This is a special case of the more general `hypot(u,v) := sqrt(u^2 + v^2)`
##'  (more details)
##' @title Compute sqrt(1 + z^2) in a numerical stable way,
##' @param z numeric or complex (or "mpfr" ..) vector (or matrix, array ..)
##' @return Numeric/complex/... vector (or array ..) with the same attributes as \code{z}
##' @author Martin Maechler
sqrt1pSqr <- function(z) {
    if(!length(z)) return(z)
    z2 <- z^2
    u <- 1+z2
    r <- sqrt(u) # direct form for normal case
    if(length(sml <- which(1 - z2^2/8 == 1))) { # also "works" for mpfr
        z22 <- z2[sml]/2
        r[sml] <- 1 + z22*(1 - z22/2) # 2-term approx 1 + z^2/2 - z^4/8
    }
    r
}
## TODO: check etc !!!! ==> ~/R/MM/NUMERICS/sqrt1pSqr-hypot.R
## ====                     ---------------------------------

## "FIXME": From ~/R/MM/NUMERICS/complex.R ,  15 Jan 2000 ((c) Martin Maechler):
##  -----   This will be better & faster __only__ in the  is.numeric() case
hypot <- function(x,y) Mod(x + 1i*y)


## NB: Binomial series of (1+u)^(1/2), exact coefficients  c_j / 2^(2j-1) :
if(FALSE) {
    k <- 0:14; noquote(formatC(choose(1/2, k) * 2^(2*k-1), w=1, format="fg"))
    ##  [1] 0.5     1       -1      2       -5      14      -42     132     -429
    ## [10] 1430    -4862   16796   -58786  208012  -742900
}
