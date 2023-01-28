## Get  dnt  from 'DPQ' and auto-produce
## i.e., "mpfr-ize" = mpfrize, see also ./dpq-mpfrized.R

## dntJKBm() := a version for Rmpfr -- by replacing the final "return line":
nB <- length(bdy <- body(dntJKBm <- dntJKBf))
new <- FALSE # to prevent FP in codetools checks
body(dntJKBm)[[nB]] <- substitute(new("mpfr", THIS), list(THIS = bdy[[nB]]))
##                                     -----------

rm(nB, bdy, new)

