require(DPQmpfr)
require(Rmpfr)

## very accurate set of *large* abscissa values
x <- mpfr(2,1024)^ seq(3, 100, by=1/2)
xN <- asNumeric(x)

## Looking at pnormU_S53() = upper bound to pnorm(), for large x,
## the term  log(4 / (3 + sqrt(1 + 8/x^2))) early underflows to 0:
eps_S53  <- function(x, method = c("direct", "log1pExct", "deltaFull",
                                   "delta1", "delta2", "delta3")) {
    method <- match.arg(method)
    delta <- function(x) { xm2 <- x^-2; (.5*(1 - sqrt(1 + 8*xm2)) + xm2) / (1 + xm2) }
    switch(method,
           "direct"    = log(4 / (3 + sqrt(1 + 8*x^-2))),
           "log1pExct" = log1p(delta(x)),
           "deltaFull" = delta(x),
           "delta1"    = -1/(x^2+1) ,               ## 1 term  for delta(x)
           "delta2"    = -(1 - (2/x)^2) / (x^2+1) , ## 2 terms for delta(x)
           "delta3"    = {tx2 <- (2/x)^2 ;  -(1 - tx2*(1 - tx2)) / (x^2+1) }, ## 3 terms
           stop("invalid 'method': ", method))
}


(S53meth <- eval(formals(eps_S53)$method))
##  "direct" "log1pExct"  "deltaFull"  "delta1" "delta2" "delta3"

epsAllS53 <- function(x, meth = eval(formals(eps_S53)$method)) {
    stopifnot(is.character(meth), is.numeric(x) || is(x,"mpfr"))
    Rl <- lapply(meth, function(m) eps_S53(x, method=m))
    new("mpfrMatrix", unlist(Rl),
        Dim = c(length(x), length(meth)),
        Dimnames = list(NULL, meth))
}

EA <- epsAllS53(x)
asNumeric(EA) # looks +- ok

## testing cbind() i.e., the Rmpfr-method  a bit here !
chkDn <- function(r) stopifnot(identical(
                         dimnames(r),
                         list(NULL, c("x", "log2x", eval(formals(eps_S53)$method)))))
chkDn2 <- function(r) stopifnot(identical(
                          colnames(r)[-(1:2)], eval(formals(eps_S53)$method)))
chkDn ( r1 <- cbind(x, log2x = log2(x), EA))
chkDn2( r2 <- selectMethod("cbind","mpfr")(x, log2x = log2(x), EA) )
## are the same if names are *made* to match:
stopifnot(identical(r1, local({t <- r2; dimnames(t) <- dimnames(r1); t})))
chkDn2( r3 <- selectMethod("cbind","mpfr")(x, log2(x), epsAllS53(x)) )
chkDn2( r4 <- cbind(x, log2(x), epsAllS53(x)) )
stopifnot(identical(r1, local({t <- r4; dimnames(t) <- dimnames(r1); t})))

U_S53mat <- function(x) cbind(x, log2x = log2(x), epsAllS53(x))
chkDn(Rx <- U_S53mat(x))
dimnames(Rx) # .. all empty, i.e.,  list(NULL, NULL)  ... even though epsAll..() gives nice colnames
##
if(FALSE) ## no longer needed !
## A version with correct colnames :
U_S53mat <- function(x) {
    eA <- epsAllS53(x)
    r <- cbind(x, log2(x), eA)
    colnames(r) <- c("x", "log2.x", colnames(eA))
    r
}

stopifnot(Rx[,-(1:2)] < 0)

op <- options(width = 100)
print(asNumeric(Rx), digits=5) ; options(op)

matplot(Rx[,"x"], - Rx[,-(1:2)], type = "b", log="xy")

## not useful really
opa <- par(ask=TRUE)
for(j in 3:ncol(Rx)) {
    D <- Rx[,-c(1:2,j)] - Rx[,j]
    neg <- mean(D < 0) > 1/2
    matplot(Rx[,"x"], if(neg) -D else D,
            main = sprintf("difference %s(Rx[,*] - Rx[, \"%s\"])",
                           if(neg) " -" else "", colnames(Rx)[j]),
            type = "b", log="xy")
}
par(opa)
## rather:
str(epsMn <- rowMeans(Rx[,-(1:2)]))

matplot(Rx[,"x"],     Rx[,-(1:2)] - epsMn, type = "b", log="x")
matplot(Rx[,"x"], abs(Rx[,-(1:2)] - epsMn), type = "b", log="xy") # don't see really

print(asNumeric(cbind(Rx[, 1:2], Rx[,-(1:2)] - epsMn)), digits=4)




if(FALSE) ## not used here anymore :
Mx <- cbind(x, log2x = log2(x), hSqrt = .5*(1-sqrt(1 + 8/x^2)),  T1 = -2*x^-2, D = -2*x^-2 - .5*(1-sqrt(1 + 8/x^2)), D2= -2*x^-2*(1 - 2*x^-2) - .5*(1-sqrt(1 + 8/x^2)))



