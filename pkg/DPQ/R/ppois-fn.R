ppoisD <- function(q, lambda)
{
  ## Purpose: ppois() via direct computation ("sum(dpois)"
  ## ----------------------------------------------------------------------
  ## Arguments: as ppois(); not yet "lower.tail" or "log.p"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  2 Mar 2004; .C -> .Call(): 2019-08-15

  .Call(C_ppoisD, as.double(q), lambda)
}


ppoisErr <- function(lambda, ppFUN = ppoisD, iP = 1e-15, verbose = FALSE) {
    ## Purpose: Given lambda, find "worst case 'x' for ppois(x, lambda)
    ##		and return the relative error of ppois()
    ## ----------------------------------------------------------------------
    ## Arguments: lambda:  parameter lambda
    ##		  iP : 1 - p; where maximal x is the p-th quantile;
    ##			as we find, the worst x we need is for p ~= 0.229
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 1 Mar 2004;  Aug. 2019 ('ppFUN' as argument)

    stopifnot(is.function(ppFUN),
              length(ff <- formals(ppFUN)) >= 2, names(ff)[[2]] == "lambda")
    ## NB: x will be large for large lambda! (R >= 3.5.0: uses 'compact' ALTREP)
    x <- 0:qpois(1-iP, lambda = lambda)
    i0 <- which.max(abs((p1x <- ppois(x, lambda = lambda)) -
                        (p0x <- ppFUN(x, lambda = lambda))))
    x0 <-   x[i0] # = i0 - 1
    p1 <- p1x[i0] # = ppois(x0, lambda)
    p2 <- p0x[i0] # ~= sum(dpois(x0:0, lambda)) # supposedly "true" value
    if(verbose)
        cat("x0=",formatC(x0),", ppois()=",formatC(p1),"\n")
    ## return rel.error *and* 'x0' as attribute :
    structure(1 - p1/p2, x0=x0)
}
