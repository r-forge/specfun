\name{p1l1}
\title{Numerically Stable p1l1(t) = (t+1)*log(1+t) - t}
\alias{p1l1}
\alias{p1l1.}
\alias{p1l1p}
\alias{p1l1ser}
\description{
  The binomial deviance function \code{\link{bd0}(x,M)} can mathematically
  be re-written as \eqn{bd0(x,M) = M * p1l1((x-M)/M)} where we look into providing
  numerically stable formula for \eqn{p1l1(t)} as it's mathematical formula
  \eqn{p1l1(t) = (t+1)\log(1+t) - t}{p1l1(t) = (t+1)*log(1+t) - t}
  suffers from cancellation for small \eqn{|t|}, even when
  \code{\link{log1p}(t)} is used instead of \code{log(1+t)}.

  Using a hybrid implementation, \code{p1l1()} uses a direct formula, now
  the stable one in \code{p1l1p()}, for \eqn{\left| t \right| > c}{|t| > c}
  and a series approximation for \eqn{\left|t\right| \le c}{|t| <= c} for
  some \eqn{c}.

  NB:  The re-expression \code{\link{log1pmx}()} is almost perfect; it
  fixes the cancellation problem entirely (and exposes the fact that
  \code{\link{log1pmx}()}'s internal cutoff seems sub optimal.
}
\usage{% -->  ../R/dgamma.R <<<<<
p1l1p  (t, ...)
p1l1.  (t)
p1l1   (t,    F = t^2/2)
p1l1ser(t, k, F = t^2/2)
}
\arguments{
  \item{t}{numeric vector ("mpfr" included), larger (or equal) to -1.}
  \item{...}{optional (tuning) arguments, passed to \code{\link{log1pmx}()}.}
  \item{k}{small positive integer, the number of terms to use in the Taylor
    series approximation \code{p1l1ser(t,k)} of \code{p1l1(t)}.}
  \item{F}{numeric vector of multiplication \bold{f}actor; \emph{must} be
    \code{t^2/2} for the \code{p1l1()} function, but can be modified,
    e.g. in more direct \code{bd0()} computations.}
}
\details{
  for now see in \code{\link{bd0}()}.
}
\value{
  numeric vector \dQuote{as} \code{x}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{bd0}};
  \code{\link{dbinom}} the latter for the C.Loader(2000) reference.
}
\examples{
t <- seq(-1, 4, by=1/64)
plot(t, p1l1ser(t, 1), type="l")
lines(t, p1l1.(t), lwd=5, col=adjustcolor(1, 1/2)) # direct formula
for(k in 2:6) lines(t, p1l1ser(t, k), col=k)

## zoom in
t <- 2^seq(-59,-1, by=1/4)
t <- c(-rev(t), 0, t)
stopifnot(!is.unsorted(t))
k.s <- 1:12; names(k.s) <- paste0("k=", 1:12)

## True function values: use Rmpfr with 256 bits precision: ---
### eventually move this to ../tests/ & ../vignettes/log1pmx-etc.Rnw
#### FIXME: eventually replace with  if(requireNamespace("Rmpfr")){ ......}
#### =====
if((needRmpfr <- is.na(match("Rmpfr", (srch0 <- search())))))
    require("Rmpfr")
p1l1.T <- p1l1.(mpfr(t, 256)) # "true" values
p1l1.n <- asNumeric(p1l1.T)
p1tab <-
    cbind(b1 = bd0(t+1, 1),
          b.10 = bd0(10*t+10,10)/10,
          dirct = p1l1.(t),
          p1l1p = p1l1p(t),
          p1l1  = p1l1 (t),
          sapply(k.s, function(k) p1l1ser(t,k)))
matplot(t, p1tab, type="l", ylab = "p1l1*(t)")
## (absolute) error:
##' legend for matplot()
mpLeg <- function(leg = colnames(p1tab), xy = "top", col=1:6, lty=1:5, lwd=1,
                  pch = c(1L:9L, 0L, letters, LETTERS)[seq_along(leg)], ...)
    legend(xy, legend=leg, col=col, lty=lty, lwd=lwd, pch=pch, ncol=3, ...)

titAbs <- "Absolute errors of p1l1(t) approximations"
matplot(t, asNumeric(p1tab - p1l1.T), type="o", main=titAbs); mpLeg()
i <- abs(t) <= 1/10 ## zoom in a bit
matplot(t[i], abs(asNumeric((p1tab - p1l1.T)[i,])), type="o", log="y",
        main=titAbs, ylim = c(1e-18, 0.003)); mpLeg()
## Relative Error
titR <- "|Relative error| of p1l1(t) approximations"
matplot(t[i], abs(asNumeric((p1tab/p1l1.T - 1)[i,])), type="o", log="y",
        ylim = c(1e-18, 2^-10), main=titR)
mpLeg(xy="topright", bg= adjustcolor("gray80", 4/5))
i <- abs(t) <= 2^-10 # zoom in more
matplot(t[i], abs(asNumeric((p1tab/p1l1.T - 1)[i,])), type="o", log="y",
        ylim = c(1e-18, 1e-9))
mpLeg(xy="topright", bg= adjustcolor("gray80", 4/5))


## Correct number of digits
corDig <- asNumeric(-log10(abs(p1tab/p1l1.T - 1)))
cbind(t, round(corDig, 1))# correct number of digits

matplot(t, corDig, type="o", ylim = c(1,17))
(cN <- colnames(corDig))
legend(-.5, 14, cN, col=1:6, lty=1:5, pch = c(1L:9L, 0L, letters), ncol=2)

## plot() function >>>> using global (t, corDig) <<<<<<<<<
p.relEr <- function(i, ylim = c(11,17), type = "o",
                    leg.pos = "left", inset=1/128,
                    main = sprintf(
                        "Correct #{Digits} in p1l1() approx., notably Taylor(k=1 .. \%d)",
                                   max(k.s)))
{
    if((neg <- all(t[i] < 0)))
        t  <- -t
    stopifnot(all(t[i] > 0), length(ylim) == 2) # as we use log="x"
    matplot(t[i], corDig[i,], type=type, ylim=ylim, log="x", xlab = quote(t), xaxt="n",
            main=main)
    legend(leg.pos, cN, col=1:6, lty=1:5, pch = c(1L:9L, 0L, letters), ncol=2,
           bg=adjustcolor("gray90", 7/8), inset=inset)
    t.epsC <- -log10(c(1,2,4)* .Machine$double.eps)
    axis(2, at=t.epsC, labels = expression(epsilon[C], 2*epsilon[C], 4*epsilon[C]),
         las=2, col=2, line=1)
    tenRs <- function(t) floor(log10(min(t))) : ceiling(log10(max(t)))
    tenE <- tenRs(t[i])
    tE <- 10^tenE
    abline (h = t.epsC,
            v = tE, lty=3, col=adjustcolor("gray",.8), lwd=2)
    AX <- if(requireNamespace("sfsmisc")) sfsmisc::eaxis else axis
    AX(1, at= tE, labels = as.expression(
                      lapply(tenE,
                             if(neg)
                                 function(e) substitute(-10^{E}, list(E = e+0))
                             else
                                 function(e) substitute( 10^{E}, list(E = e+0)))))
}

p.relEr(t > 0, ylim = c(1,17))
p.relEr(t > 0) # full positive range
p.relEr(t < 0) # full negative range
if(FALSE) {## (actually less informative):
 p.relEr(i = 0 < t & t < .01)  ## positive small t
 p.relEr(i = -.1 < t & t < 0) ## negative small t
}

## Find approximate formulas for accuracy of k=k*  approximation
d.corrD <- cbind(t=t, as.data.frame(corDig))
names(d.corrD) <- sub("k=", "nC_",  names(d.corrD))

fmod <- function(k, data, cut.y.at = -log10(2 * .Machine$double.eps),
                 good.y = -log10(.Machine$double.eps), # ~ 15.654
                 verbose=FALSE) {
    varNm <- paste0("nC_",k)
    stopifnot(is.numeric(y <- get(varNm, data, inherits=FALSE)),
              is.numeric(t <- data$t))# '$' works for data.frame, list, environment
    i <- 3 <= y & y <= cut.y.at
    i.pos <- i & t > 0
    i.neg <- i & t < 0
    if(verbose) cat(sprintf("k=\%d >> y <= \%g ==> #{pos. t} = \%d ;  #{neg. t} = \%d\n",
                            k, cut.y.at, sum(i.pos), sum(i.neg)))
    nCoefLm <- function(x,y) `names<-`(.lm.fit(x=x, y=y)$coeff, c("int", "slp"))
    nC.t <- function(x,y) { cf <- nCoefLm(x,y); c(cf, t.0 = exp((good.y - cf[[1]])/cf[[2]])) }
    cbind(pos = nC.t(cbind(1, log( t[i.pos])), y[i.pos]),
          neg = nC.t(cbind(1, log(-t[i.neg])), y[i.neg]))
}
rr <- sapply(k.s, fmod, data=d.corrD, verbose=TRUE, simplify="array")
stopifnot(rr["slp",,] < 0) # all slopes are negative (important!)
matplot(k.s, t(rr["slp",,]), type="o", xlab = quote(k), ylab = quote(slope[k]))
## fantastcally close to linear in k
## The numbers, nicely arranged
ftable(aperm(rr, c(3,2,1)))
signif(t(rr["t.0",,]),3) # ==> Should be boundaries for the hybrid p1l1()
##           pos      neg
## k=1  6.60e-16 6.69e-16
## k=2  3.65e-08 3.65e-08
## k=3  1.30e-05 1.32e-05
## k=4  2.39e-04 2.42e-04
## k=5  1.35e-03 1.38e-03
## k=6  4.27e-03 4.34e-03
## k=7  9.60e-03 9.78e-03
## k=8  1.78e-02 1.80e-02
## k=9  2.85e-02 2.85e-02
## k=10 4.13e-02 4.14e-02
## k=11 5.62e-02 5.64e-02
## k=12 7.24e-02 7.18e-02

###------------- Well,  p1l1p()  is really basically good enough ... with a small exception:
rErr1k <- curve(asNumeric(p1l1p(x) / p1l1.(mpfr(x, 4096)) - 1), -.999, .999,
                n = 4000, col=2, lwd=2)
abline(h = c(-8,-4,-2:2,4,8)* 2^-52, lty=2, col=adjustcolor("gray20", 1/4))
## well, have a "spike" at around -0.8 -- why?

plot(abs(y) ~ x, data = rErr1k, ylim = c(4e-17, max(abs(y))),
     ylab=quote(abs(hat(p)/p - 1)),
     main = "p1l1p(x) -- Relative Error wrt mpfr(*. 4096) [log]",
     col=2, lwd=1.5, type = "b", cex=1/2, log="y", yaxt="n")
sfsmisc::eaxis(2)
eps124 <-  c(1, 2,4,8)* 2^-52
abline(h = eps124, lwd=c(3,1,1,1), lty=c(1,2,2,2), col=adjustcolor("gray20", 1/4))
axLab <- expression(epsilon[c], 2*epsilon[c], 4*epsilon[c], 8*epsilon[c])
axis(4, at = eps124, labels = axLab, col="gray20", las=1)
abline(v= -.791, lty=3, lwd=2, col="blue4") # -.789  from visual ..
##--> The "error" is in log1pmx() which has cutoff minLog1Value = -0.79149064
##--> which is clearly not optimal, at least not for computing p1l1p()

d <- 1/2048; x <- seq(-1+d, 1, by=d)
p1l1Xct <- p1l1.(mpfr(x, 4096))
rEx.5 <- asNumeric(p1l1p(x, minL1 = -0.5) / p1l1Xct - 1)
lines(x, abs(rEx.5), lwd=2.5, col=adjustcolor(4, 1/2)); abline(v=-.5, lty=2,col=4)
rEx.25 <- asNumeric(p1l1p(x, minL1 = -0.25) / p1l1Xct - 1)
lines(x, abs(rEx.25), lwd=3.5, col=adjustcolor(6, 1/2)); abline(v=-.25, lty=2,col=6)
lines(lowess(x, abs(rEx.5),  f=1/20), col=adjustcolor(4,offset=rep(1,4)/3), lwd=3)

lines(lowess(x, abs(rEx.25), f=1/20), col=adjustcolor(6,offset=rep(1,4)/3), lwd=
3)

rEx.4 <- asNumeric(p1l1p(x, tol_logcf=1e-15, minL1 = -0.4) / p1l1Xct - 1)
lines(x, abs(rEx.4), lwd=5.5, col=adjustcolor("brown", 1/2)); abline(v=-.25, lty=2,col="brown")

if(needRmpfr && isNamespaceLoaded("Rmpfr"))
    detach("package:Rmpfr")
}% examples
