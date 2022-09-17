#### Extreme-tail pnorm() and qnorm() approximations
#### ===============================================
### notably  pnormAsymp()  and  qnormAsymp()
###          ~~~~~~~~~~~~       ~~~~~~~~~~~~
## partly inspired by MM's NUMERICS/dpq-functions/qnorm-asymptotic.R
## and also ../man/pnormAsymp.Rd
##          ../man/qnormAsymp.Rd

library(DPQ)
library(sfsmisc)# {we import it in DPQ}

r <- sort(c(seq(1,100, by=1/8),
            2 ^ c(seq(6.5, 10, by=1/16), seq(10.25, 14.5, by=1/4))))
str(r)
summary(r)
## Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
## 1.00    28.09    55.19   231.97    82.28 23170.47

s <- r^2
lp <- -s

## Start with quantiles so we know the  "truth according to pnorm()" :
qs <- 2^seq( 0, 35, by=1/256) # => s >= 1.84  --> no NaN in xs_5 etc
## == the "true" qnorm() according to pnorm()
lp  <- pnorm(qs, lower.tail=FALSE, log.p=TRUE)
qnp <- qnorm(lp, lower.tail=FALSE, log.p=TRUE)
s <- -lp # =  -log(1 - Phi(qs))
summary(r <- sqrt(s))
## store in 'dat' including current "inaccuracy" of qnorm():
head(dat <- data.frame(l2q=log2(qs), qs, # lp,
                       s,  # r = sqrt(s), t = log(s),
                       relE_qn = relErrV(qs, qnp)))
				##^^^^^^  will depend much on R version

p.ver <- function() mtext(R.version.string, cex=3/4, adj=1)
readUser <- function(i, max.i) {
    if(interactive() && i != max.i) {
        cat("[Enter] to continue: ")
        cat(readLines(stdin(), n=1), "\n")
    }
}


## relative Error (log x)
plot(relE_qn ~ s, dat, log="x", type="l", col = 2); p.ver()

## |rel.Error| ==> log-log scale
plot(abs(relE_qn) ~ s, dat, log="xy", type="l", col = 2, axes=FALSE)
eaxis(1); eaxis(2); p.ver()

## now look at the asymptotic approximations:
k.s <- 0:5; nks <- paste0("k=", k.s)
qnAsy <- sapply(setNames(k.s, nks), function(ord)
    qnormAsymp(lp, lower.tail=FALSE, log.p=TRUE, order=ord))
relEAsy <- qnAsy / qs - 1

matplot(sqrt(s), relEAsy, type="l", log="x", xlab = quote(r == sqrt(s)),
        main = "relative Error of qnormAsymp(.., k=*)"); p.ver()
legend("top", nks, col=1:6, lty=1:5, bty="n")

mp <- c(-1,1)
xL <- list(NULL, c(2, 1000),  c(2, 1e5), c(5, 1e5), c(10, 1e8), c(10, 1e10))
yL <- list(NULL,   mp * 0.1,  mp * 1e-4, mp * 1e-7, mp * 1e-11, mp * 2e-15)
stopifnot(length(xL) == length(yL))
for(j in seq_along(xL)) {
    matplot(sqrt(s), relEAsy, type="l", log="x", xaxt="n", lwd=2,
            xlab = quote(r == sqrt(s)),
            xlim = xL[[j]], ylim = yL[[j]],
            main = "relative Error of qnormAsymp(.., k=*)")
    eaxis(1, sub10=2) # not here: p.ver()
    legend("top", nks, col=1:6, lty=1:5, lwd=2, bty="n")
    readUser(j, length(xL))
}

absP <- function(re) pmax(abs(re), 2e-17) # not zero, so log-scale "shows" it

for(j in seq_along(xL)) {
    if(!is.null(yli <- yL[[j]])) yli <- c(2e-17, 100*yli[2])
    if(!is.null(xli <- xL[[j]])) xli[2] <- 10*xli[2]
    cat("j=",j,"; rbind(xlim,ylim) = \n"); print(rbind(xlim=xli, ylim=yli))
    matplot(sqrt(s), absP(relEAsy), type="l", log="xy", xaxt="n", yaxt="n", lwd=2,
            xlab = quote(r == sqrt(s)),
            xlim = xli, ylim = yli,
            main = "relative Error of qnormAsymp(.., k=*)")
    eaxis(1, sub10=2) # not here: p.ver()
    eaxis(2, sub10=c(-3,2))
    legend("top", nks, col=1:6, lty=1:5, lwd=2, bty="n")
    readUser(j, length(xL))
}
