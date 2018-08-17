ppoisD <- function(q, lambda)
{
  ## Purpose: ppois() via direct computation ("sum(dpois)"
  ## ----------------------------------------------------------------------
  ## Arguments: as ppois(); now "lower.tail" or "log.p"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  2 Mar 2004, 09:57

    if(length(lambda) != 1)
        stop("argument 'lambda' must have length 1 !")
    n <- length(x <- as.double(q))
    .C("ppois_D",
       x,
       n,
       as.double(lambda),
       pr = double(n),
     , PACKAGE = "DPQ")$pr
}


err.lambd0 <- function(lam, iP = 1e-15) {
    ## Purpose: Given lambda, find "worst case 'x' for ppois(x, lambda)
    ##		and return the relative error of ppois()
    ## ----------------------------------------------------------------------
    ## Arguments: lam:  parameter lambda
    ##		  iP : 1 - p; where maximal x is the p-th quantile;
    ##			as we find, the worst x we need is for p ~= 0.229
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  1 Mar 2004, 17:40

    x <- seq(0, qpois(1-iP, lambda = lam), by = 1) # large for large lambda!
    i0 <- which.max(abs(ppois(x, lambda = lam) -
                        cumsum(dpois(x, lambda = lam))))
    x0 <- i0 -1 ## == x[i0]
    p1 <- ppois(x0, lam)
    cat("x0=",formatC(x0),", ppois()=",formatC(p1),"\n")
    p2 <- sum(dpois(0:x0, lam))
    ## rel.error:
    1 - p1/p2
}

err.lambd1 <- function(lam, iP = 1e-15) {
    ## Purpose: Given lambda, find "worst case 'x' for ppois(x, lambda)
    ##		and return the relative error of ppois()
    ## ----------------------------------------------------------------------
    ## Arguments: lam:  parameter lambda
    ##		  iP : 1 - p; where maximal x is the p-th quantile;
    ##			as we find, the worst x we need is for p ~= 0.229
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  2 Mar 2004, 10:00

    x <- seq(0, qpois(1-iP, lambda = lam), by = 1) # large for large lambda!
    i0 <- which.max(abs(ppois (x, lambda = lam) -
                        ppoisD(x, lambda = lam)))
    ##	                ~~~~~~ much slower than cumsum(dpois(x,lambda=lam))  {why?}
    x0 <- i0 -1 ## == x[i0]
    p1 <- ppois(x0, lam)
    cat("x0=",formatC(x0),", ppois()=",formatC(p1),"\n")
    p2 <- sum(dpois(0:x0, lam))
    ## rel.error:
    1 - p1/p2
}

err.lambda <- err.lambd1 ## a tiny bit more accurate but __MUCH__ slower
err.lambda <- err.lambd0 ## for speed (and reproducability)
