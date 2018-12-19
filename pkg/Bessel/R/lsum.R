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
