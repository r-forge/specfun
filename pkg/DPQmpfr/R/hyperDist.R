##' dhyper() as "big rational"
##' NB: chooseZ(n, k)  can have bigz 'n' but must have integer 'k'
##' ==>  {m, n} may be bigz,  but {x, k} must be integer here :
dhyperQ <- function(x, m, n, k) {
  stopifnot(k-x == as.integer(k-x))
  chooseZ(m, x) * chooseZ(n, k-x) / chooseZ(m+n, k)
}

## Imported from 'DPQ': support of the hypergeometric distrib. as function of its parameters:
## .suppHyper <- function(m,n,k) max(0, k-n) : min(k, m)

##' phyper() -- simple version: length(x) == 1
phyperQ <- function(x, m, n, k, lower.tail=TRUE) {
    stopifnot(length(m) == 1, length(n) == 1, length(k) == 1,
              k-x == as.integer(k-x),
              (xi <- as.integer(x)) == x, length(xi) == 1)
    ## max(.) and min() such that the *first* argument may be bigz
    if(xi <  (x.min <- as.integer(max(k-n, 0L)))) return(if(lower.tail) 0 else 1)
    if(xi >= (x.max <- as.integer(min( m , k )))) return(if(lower.tail) 1 else 0)
    ## x.min <= x < x.max : Still got a feeling we should *cache* binomial coeff?
    sum( dhyperQ(if(lower.tail) x.min:xi else x.max:(xi+1L),
                 m, n, k) )
}

## really the cumsum(.) version is mostly much more useful:

## have imported from 'DPQ': .suppHyper <- function(m,n,k) max(0, k-n) : min(k, m)
## but here "plug in":
phyperQall <- function(m, n, k, lower.tail=TRUE) {
    stopifnot(length(m) == 1, length(n) == 1, length(k) == 1)
    ## in both cases, we sum up to 1 ("proving" consistency), but start positive:
    if(lower.tail)
                  cumsum(dhyperQ(max(0, k-n) : min(k, m), m,n,k))
    else
      rev.default(cumsum(dhyperQ(min(k, m) : max(0, k-n), m,n,k)))
}
