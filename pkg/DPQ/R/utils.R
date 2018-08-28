## Convenient argument checking

all_mpfr <- function(...) {
    for(x in list(...)) if(!inherits(x, "mpfr")) return(FALSE)
    ## else return
    TRUE
}

any_mpfr <- function(...) {
    for(x in list(...)) if(inherits(x, "mpfr")) return(TRUE)
    ## else return
    FALSE
}


### form01.prec  was in  source("~/R/MM/MISC/Util.R")
form01.prec <- function(x, digits = getOption("digits"), width = digits + 2,
                        eps = 1e-6, ...,
                        fun = function(x,...) formatC(x, flag='-',...))
{
  ## Purpose: format numbers in [0,1] with "precise" result,
  ##          using "1-.." if necessary.
  ## -------------------------------------------------------------------------
  ## Arguments: x:      numbers in [0,1]; (still works if not)
  ##            digits, width: number of digits, width to use with 'fun'
  ##            eps:    Use '1-' iff  x in  (1-eps, 1] -- 1e-6 is OPTIMAL
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14 May 97, 18:07
  if(as.integer(digits) < 4) stop('digits must be >= 4')
  if(eps < 0 || eps > .1) stop('eps must be in [0, .1]')
  i.swap <- 1-eps < x  &  x <= 1 #-- Use "1- ." if <~ 1,  normal 'fun' otherwise
  r <- character(length(x))
  if(any(i.swap))
    r[i.swap] <-
      paste("1-", fun(1-x[i.swap], digits=digits - 5, width=width-2, ...),
            sep='')# -5: '1-' + 4 for exponent -1 for '0' (in other case)
  if(any(!i.swap))
    r[!i.swap] <- fun(x[!i.swap], digits=digits, width=width,...)
  attributes(r) <- attributes(x)
  r
}

