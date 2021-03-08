#### Testing  dbinom_raw(), dnbinomR() and dnbinom.mu() >>> ../R/dbinom-nbinom.R <<<
####        ../man/dbinom_raw.Rd & ../man/dnbinomR.Rd
require(DPQ)

## "FIXME:" get  Matrix functions for relErr etc

###  dbinom() vs  dbinom.raw() :

for(n in 1:20) {
    cat("n=",n," ")
    for(x in 0:n)
        cat(".")
        for(p in c(0, .1, .5, .8, 1)) {
            stopifnot(all.equal(dbinom_raw(x, n, p, q=1-p, log=FALSE),
                                dbinom    (x, n, p,        log=FALSE)),
                      all.equal(dbinom_raw(x, n, p, q=1-p, log =TRUE),
                                dbinom    (x, n, p,        log =TRUE)))
    }
    cat("\n")
}

###  dnbinom*() :
stopifnot(exprs = {
    dnbinomR(0, 1, 1) == 1
})

### exploring 'eps' == "true" tests must be done with  Rmpfr !!
