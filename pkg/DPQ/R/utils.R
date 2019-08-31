## Not exported, used to make  'R CMD check <pkg>'  be faster *or* more extensive:
doExtras <- function(int = interactive()) {
    int || nzchar(Sys.getenv("R_DPQ_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

## Convenient argument checking

all_mpfr <- function(...) {
    ## Remain lazy ==> do *NOT* use list(...)  which evaluates all
    for(i in ...length()) if(!inherits(...elt(i), "mpfr")) return(FALSE)
    ## else return
    TRUE
}

any_mpfr <- function(...) {
    ## Remain lazy ==> do *NOT* use list(...)  which evaluates all
    for(i in ...length()) if(inherits(...elt(i), "mpfr")) return(TRUE)
    ## else return
    FALSE
}

## MM: From ~/R/D/r-devel/R/src/nmath/pgamma.c : --- Help page now  >>> ../man/logspace.add.Rd <<<<
logspace.add <- function(lx, ly) pmax(lx, ly) + log1p(exp(-abs(lx - ly)))

logspace.sub <- function(lx, ly) lx + log1mexp(lx - ly)



### was as 'form01.prec' in  source("~/R/MM/MISC/Util.R")
format01prec <- function(x, digits = getOption("digits"), width = digits + 2,
                        eps = 1e-6, ...,
                        FUN = function(x,...) formatC(x, flag='-',...))
{
  ## Purpose: format numbers in [0,1] with "precise" result,
  ##          using "1-.." if necessary.
  ## -------------------------------------------------------------------------
  ## Arguments: x:      numbers in [0,1]; (still works if not)
  ##            digits, width: number of digits, width to use with 'FUN'
  ##            eps:    Use '1-' iff  x in  (1-eps, 1] -- 1e-6 is OPTIMAL
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14 May 97, 18:07
  if(as.integer(digits) < 4) stop('digits must be >= 4')
  if(eps < 0 || eps > .1) stop('eps must be in [0, .1]')
  i.swap <- 1-eps < x  &  x <= 1 #-- Use "1- ." if <~ 1,  normal 'FUN' otherwise
  r <- character(length(x))
  if(hasNA <- anyNA(i.swap)) {
    iNA <- is.na(i.swap)
    i.swap  <-  i.swap & !iNA
    ni.swap <- !i.swap & !iNA
  } else
    ni.swap <- !i.swap
  if(any(i.swap))
      r[i.swap] <- paste0("1-",
                          FUN(1-x[i.swap],
                              ## -5: '1-' + 4 for exponent -1 for '0' (in other case)
                              digits=digits - 5, width=width-2, ...))
  if(any(ni.swap))
    r[ni.swap] <- FUN(x[ni.swap], digits=digits, width=width,...)
  if(hasNA)
    r[iNA] <- "NA"
  attributes(r) <- attributes(x)
  r
}

