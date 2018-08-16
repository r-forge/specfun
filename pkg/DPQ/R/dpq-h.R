## Define S functions for the utility macros in
### ~/R/D/r-devel/R/src/nmath/dpq.h
### ---------------------------------
## if(FALSE) ## 'include' this via
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")

.D_0 <- function(log.p) if(log.p) -Inf else 0
.D_1 <- function(log.p) if(log.p)   0	 else 1

##' == 0 when (lower.tail=TRUE, log.p=FALSE)
.DT_0 <- function(lower.tail, log.p)
    if(lower.tail) .D_0(log.p) else .D_1(log.p)

##' == 1 when (lower.tail=TRUE, log.p=FALSE)
.DT_1 <- function(lower.tail, log.p)
    if(lower.tail) .D_1(log.p) else .D_0(log.p)

.D_Lval <- function(p, lower.tail) if(lower.tail) p else (1 - p) #   p
.D_Cval <- function(p, lower.tail) if(lower.tail) (1 - p) else p # 1 - p


.D_val <- function(x, log.p)  if(log.p) log(x) else x	     # x  in pF(x,..)
.D_qIv <- function(p, log.p)  if(log.p) exp(p) else p	     # p  in qF(p,..)
.D_exp <- function(x, log.p)  if(log.p) x else exp(x)	     # exp(x)
.D_log <- function(p, log.p)  if(log.p) p else log(p)	     # log(p)
.D_Clog<- function(p, log.p)  if(log.p) log1p(-p) else ((0.5 - p) + 0.5)# [log](1-p)

M.LN2 <- log(2)

## log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x))  == log1mexp(-x)
R.Log1.Exp <- function(x) ifelse(x > -M.LN2, log(-expm1(x)), log1p(-exp(x)))

## log(1-exp(x)): R_D_LExp(x) == (log1p(- .D_qIv(x))) but even more stable:
.D_LExp <- function(x, log.p) if(log.p) R.Log1.Exp(x) else log1p(-x)

.DT_val  <- function(x, lower.tail, log.p)
    .D_val(.D_Lval(x, lower.tail), log.p) #	x   in pF
.DT_Cval <- function(x, lower.tail, log.p)
    .D_val(.D_Cval(x, lower.tail), log.p) # 1 - x in pF

## .DT_qIv <- function(p)	.D_Lval(.D_qIv(p))	  #   p	 in qF !
.DT_qIv <- function(p, lower.tail, log.p) {
    if(log.p) if(lower.tail) exp(p) else - expm1(p)
    else .D_Lval(p, lower.tail)
}

## .DT_CIv <- function(p)	.D_Cval(.D_qIv(p))	  #  1 - p in qF
.DT_CIv <- function(p, lower.tail, log.p) {
    if(log.p) if(lower.tail) -expm1(p) else exp(p)
    else .D_Cval(p, lower.tail)
}

.DT_exp  <- function(x, lower.tail, log.p)		# exp( x )
    .D_exp(.D_Lval(x, lower.tail), log.p)

.DT_Cexp <- function(x, lower.tail, log.p)		# exp(1 - x)
    .D_exp(.D_Cval(x, lower.tail), log.p)

.DT_log <- function(p, lower.tail, log.p)		# log (p ) in qF
    if(lower.tail) .D_log(p, log.p) else .D_LExp(p, log.p)

.DT_Clog <- function(p, lower.tail, log.p)		# log(1-p) in qF
    if(lower.tail) .D_LExp(p, log.p) else .D_log(p, log.p)

.DT_Log <- function(p, lower.tail)		# == .DT_log when we "know" log.p=TRUE
    if(lower.tail) p else R.Log1.Exp(p)		# log(p) in qF	[for log_p = TRUE ]


##  R_Q_P01_boundaries <- function(p, _LEFT_, _RIGHT_)
##  ------------------
##  cannot work: in C they are macros where return(.) is toplevel !