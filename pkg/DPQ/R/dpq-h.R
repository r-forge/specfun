## Define S functions for the utility macros in
### ~/R/D/r-devel/R/src/nmath/dpq.h
### ---------------------------------
## if(FALSE) ## 'include' this via
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")

R.D..0 <- function(log.p) if(log.p) -Inf else 0
R.D..1 <- function(log.p) if(log.p)   0	 else 1

##' == 0 when (lower.tail=TRUE, log.p=FALSE)
R.DT.0 <- function(lower.tail, log.p)
    if(lower.tail) R.D..0(log.p) else R.D..1(log.p)

##' == 1 when (lower.tail=TRUE, log.p=FALSE)
R.DT.1 <- function(lower.tail, log.p)
    if(lower.tail) R.D..1(log.p) else R.D..0(log.p)

R.D.Lval <- function(p, lower.tail) if(lower.tail) p else (1 - p) #   p
R.D.Cval <- function(p, lower.tail) if(lower.tail) (1 - p) else p # 1 - p


R.D.val <- function(x, log.p)  if(log.p) log(x) else x	     # x  in pF(x,..)
R.D.qIv <- function(p, log.p)  if(log.p) exp(p) else p	     # p  in qF(p,..)
R.D.exp <- function(x, log.p)  if(log.p) x else exp(x)	     # exp(x)
R.D.log <- function(p, log.p)  if(log.p) p else log(p)	     # log(p)
R.D.Clog<- function(p, log.p)  if(log.p) log1p(-p) else ((0.5 - p) + 0.5)# [log](1-p)

M.LN2 <- log(2)

## log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x))  == log1mexp(-x)
R.Log1.Exp <- function(x) ifelse(x > -M.LN2, log(-expm1(x)), log1p(-exp(x)))

## log(1-exp(x)): R_D_LExp(x) == (log1p(- R.D.qIv(x))) but even more stable:
R.D.LExp <- function(x, log.p) if(log.p) R.Log1.Exp(x) else log1p(-x)

R.DT.val  <- function(x, lower.tail, log.p)
    R.D.val(R.D.Lval(x, lower.tail), log.p) #	x   in pF
R.DT.Cval <- function(x, lower.tail, log.p)
    R.D.val(R.D.Cval(x, lower.tail), log.p) # 1 - x in pF

## R.DT.qIv <- function(p)	R.D.Lval(R.D.qIv(p))	  #   p	 in qF !
R.DT.qIv <- function(p, lower.tail, log.p) {
    if(log.p) if(lower.tail) exp(p) else - expm1(p)
    else R.D.Lval(p, lower.tail)
}

## R.DT.CIv <- function(p)	R.D.Cval(R.D.qIv(p))	  #  1 - p in qF
R.DT.CIv <- function(p, lower.tail, log.p) {
    if(log.p) if(lower.tail) -expm1(p) else exp(p)
    else R.D.Cval(p, lower.tail)
}

R.DT.exp  <- function(x, lower.tail, log.p)		# exp( x )
    R.D.exp(R.D.Lval(x, lower.tail), log.p)

R.DT.Cexp <- function(x, lower.tail, log.p)		# exp(1 - x)
    R.D.exp(R.D.Cval(x, lower.tail), log.p)

R.DT.log <- function(p, lower.tail, log.p)		# log (p ) in qF
    if(lower.tail) R.D.log(p, log.p) else R.D.LExp(p, log.p)

R.DT.Clog <- function(p, lower.tail, log.p)		# log(1-p) in qF
    if(lower.tail) R.D.LExp(p, log.p) else R.D.log(p, log.p)

R.DT.Log <- function(p, lower.tail)		# == R.DT.log when we "know" log.p=TRUE
    if(lower.tail) p else R.Log1.Exp(p)		# log(p) in qF	[for log_p = TRUE ]


##  R_Q_P01_boundaries <- function(p, _LEFT_, _RIGHT_)
##  ------------------
##  cannot work: in C they are macros where return(.) is toplevel !
