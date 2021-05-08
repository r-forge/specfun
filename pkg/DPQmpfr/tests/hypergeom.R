library(DPQ)
library(DPQmpfr)

## Imported from 'DPQ', support of the hypergeometric distrib. as function of its parameters:
## .suppHyper <- DPQ::.suppHyper
.suppHyper <- function(m,n,k) max(0, k-n) : min(k, m)

## case where phyper(*, log=TRUE) fails :
phyper(11, 15, 0, 12) # 0 -- correct
phyper(11, 15, 0, 12, log=TRUE) # NaN  || correct value:  log(0) = -Inf
(phyp5.0.12 <- cumsum(dhyperQ(0:12, m=15,n=0,k=12)))

require(Rmpfr)
.N <- asNumeric # short form for here

m <- 1000; n <- 1001; k <- 1500
dhyp.1k <- dhyperQ(.suppHyper(m,n,k), m,n,k)
phyp.1k <- cumsum(dhyp.1k)
## Rel.Error: number of correct digits
.N(-log10(abs(1 - phyper(.suppHyper(m,n,k), m,n,k) / .bigq2mpfr(phyp.1k, 99))))

## even a bit larger
m <- 9999; n <- 500; k <- 10100
range(x <- .suppHyper(m,n,k)) # 9600 9999
dQ <- dhyperQ(x, m,n,k) # these add exactly to 1 :
stopifnot(sum(dQ) == 1)


dM <- .bigq2mpfr(dQ, precB=128)
tail( log(dM) ) # -1395.7599 ... -1443.7188576

1 - phyper(x, m,n,k, log.p=TRUE) / log(cumsum(dM)) # from 0 --> 1  ???
## -->
## this is *bad*: cumsum(dM) suffers from cancellation still!
cbind(x, phyper = phyper(x, m,n,k), cum.dQ = asNumeric(cumsum(dM)),
      phyp.log  = phyper(x, m,n,k, log.p=TRUE), log.cumQ = asNumeric(log(cumsum(dM))))

## this shows the cancellation:
cbind(x,
      cumD = .N(1 - cumsum(dM)),
      revCS= .N(rev(cumsum(rev(dM)))))[ 1:99, ]
## Using bigrational ("Q") is still exact:   ---- well there's problem somewhere !
cbind(x,
      cumD = .N(1 - cumsum(    dQ)),
      revCS= .N(rev(cumsum(rev(dQ)))))[ 1:99, ]

## seems good
I_CSd <- .bigq2mpfr(1 - cumsum(dQ), precB=1024)
roundMpfr(I_CSd, precBits = 8)
.N(log(I_CSd)) # should be "exact"

## rel.err :
relE.up.ln <- .N(1- phyper(x, m,n,k, lower.tail = FALSE, log.p = TRUE) / log(I_CSd))
plot(x, relE.up.ln, type = "o", cex = 1/2)
abline(h = c(-1,1)*.Machine$double.eps, col = "thistle", lty=2)
plot(x, log(I_CSd), type="l")
## rel.err of phyper(*, log=TRUE) : --> is exploding!
1- phyper(x, m,n,k, log.p=TRUE) / log(cumsum(dM))


