#### qpois(), qbinom() and qnbinom() overflow much too early in the upper tail
#### ~~~~~    ~~~~~~       ~~~~~~~

if(!dev.interactive(orNone=TRUE)) pdf("qPoisBinom-ex.pdf")
.O.P. <- par(no.readonly=TRUE)

## It's source code  ~/R/D/r-devel/R/src/nmath/qpois.c  contains
'
    /* Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
     * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if(!lower_tail || log_p) {
	p = R_DT_qIv(p); /* need check again (cancellation!): */
	if (p == R_DT_0) return 0;
	if (p == R_DT_1) return ML_POSINF;
    }
    /* temporary hack --- FIXME --- */
    if (p + 1.01*DBL_EPSILON >= 1.) return ML_POSINF;
'
## ==> same "bug" also in qbinom() and qnbinom()  [! ?? ]

e <- c(-1000, -500, -200, -100, seq(-70,-1/4, by=1/4))

lambda <- 10000
qp <- qpois(2^e, lambda=lambda, lower.tail=FALSE)

## 'cbind_no_rownames' :
cbNoRN <- function(...) {
    r <- cbind(...)
    dimnames(r)[[1]] <- rep.int("", nrow(r))
    r
}

cbNoRN(e, p=2^e, qpois=qp)
##         e             p qpois
##  -1000.00 9.332636e-302   Inf
##     .....  ............   ...
##     .....  ............   ...
##    -52.25  1.867165e-16   Inf
##    -52.00  2.220446e-16   Inf
##    -51.75  2.640570e-16   Inf
##    -51.50  3.140185e-16 10770
##    -51.25  3.734330e-16 10769
##    -51.00  4.440892e-16 10769
##    -50.75  5.281140e-16 10769
##    -50.50  6.280370e-16 10769
##    -50.25  7.468660e-16 10769
##    -50.00  8.881784e-16 10769

plot(qp ~ e, type = "b", subset = -(1:5),
     main = paste0("qpois(2^e, lambda=",lambda,
                   ") - early overflow to +Inf"))

## qbinom() "same" problem -- does not overflow to +Inf but to
##  'size' ( = n in C ) = 100  here : as indeed, the source
##  ~/R/D/r-devel/R/src/nmath/qbinom.c has an *ADDITIONAL* hack here :
'
    q = 1 - pr;
    if(q == 0.) return n; /* covers the full range of the distribution */
'
qBin <- qbinom(2^e, size = 100, prob = 0.4, lower.tail=FALSE)
cbNoRN(e, p=2^e, qBin)[c(1, 75:85),]

plot(qBin ~ e, type = "b", subset = -(1:5),
     main = paste0("qbinom(2^e, size = 100, prob = 0.4, lower.tail=FALSE",
                   ") - early overflow to 'size'")); abline(h=100, lty=3)

## qnbinom() "same" problem --
##  ~/R/D/r-devel/R/src/nmath/qnbinom.c
qNB <- qnbinom(2^e, size = 100, prob = 0.4, lower.tail=FALSE)
cbNoRN(e, p=2^e, qNB)[c(1, 70:82),]
 ##        e             p qNB
 ## -1000.00 9.332636e-302   0
 ##   -53.75  6.601426e-17   0
 ##   -53.50  7.850462e-17   0 <<<< !! even more wrong
 ##   -53.25  9.335826e-17 Inf
 ##   -53.00  1.110223e-16 Inf
 ##   -52.75  1.320285e-16 Inf
 ##   -52.50  1.570092e-16 Inf
 ##   -52.25  1.867165e-16 Inf
 ##   -52.00  2.220446e-16 Inf
 ##   -51.75  2.640570e-16 Inf <<
 ##   -51.50  3.140185e-16 337
 ##   -51.25  3.734330e-16 337
 ##   -51.00  4.440892e-16 337
 ##   -50.75  5.281140e-16 337

## to make the jump to Inf, then 0, more visible, replace Inf by HUGE :
qN. <- qNB
qN.[qNB == Inf] <- 1e300
plot(qN. ~ e, type = "l", subset = -(1:5), ylim = range(qNB, finite=TRUE),
     main = paste0("qnbinom(2^e, size = 100, prob = 0.4, lower.tail=FALSE)",
                   ") - early \"overflow\" to **WRONG**")); abline(h=0, lty=3)

require(DPQ)
## Embarrassing forgotten FIXME (for boring boundary cases only):
M <- 2^31; pr <- 1e-9
stopifnot(exprs = {
    qbinomR (0:1, size=M, prob=pr) == c(0, M) # was 0 1
    qnbinomR(0:1, size=M, prob=pr) == c(0, Inf) # " "
    qpoisR  (0:1,      M)          == c(0, Inf)
    qbinomR (c(-Inf,0), size=M, prob=pr, log.p=TRUE) == c(0, M)
    qnbinomR(c(-Inf,0), size=M, prob=pr, log.p=TRUE) == c(0, Inf)
    qpoisR  (c(-Inf,0),      M,          log.p=TRUE) == c(0, Inf)
})


