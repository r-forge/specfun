
###==> see also ./pnt-ex.R
###		~~~~~~~~~~
library(dcdflib)

###--- Originally from ~/R/MM/NUMERICS/dpq-functions/pnt-ex.R -----------------

stopifnot(exprs = {
    require(graphics)
    require(sfsmisc)
})

if(!dev.interactive(orNone=TRUE)) pdf("ptnc_tsts.pdf")

## This shows unpleasant behavior of ptnc(), i.e. dcdflib's cumtnc() ,
## for `small' x  --- already in version 1.1 (of dcdflib.c & cdflib.h) --
par(mfrow=c(3,2),lab=c(10,20,7))
{
    plot(function(t)ptnc(t, df = 10, ncp=5.0),0,3)# no log
    plot(function(t) pt (t, df = 10, ncp=5.0),0,3, add=TRUE, col="purple")
    abline(h=10^(-2:0),col="green",lty=3,lwd=.1)

    plot(function(t)ptnc(t, df = 10, ncp=5.0),0,3, log="y")
    ## bug: jumps to 0 for  t <= 1
    plot(function(t) pt (t, df = 10, ncp=5.0),0,3, add=TRUE, col="purple")
    abline(h=10^(-2:0),col="green",lty=3,lwd=.1)

    plot(function(t)ptnc(t, df = 10, ncp=25),0,30)# no log
    plot(function(t) pt(t, df = 10, ncp=25),0,30, add=TRUE, col="purple")
    abline(h=10^(-2:0),col="green",lty=3,lwd=.1)

    plot(function(t)ptnc(t, df = 10, ncp=25),0,30, log="y")
    ## "bug": jumps to 0 ..
    plot(function(t) pt(t, df = 10, ncp=25),0,30, add=TRUE, col="purple")
    abline(h=10^(-2:0),col="green",lty=3,lwd=.1)

    plot(function(t)ptnc(t, df = 10, ncp=35),0,70)# no log
    plot(function(t) pt(t, df = 10, ncp=35),0,70, add=TRUE, col="purple")
    abline(h=10^(-2:0),col="green",lty=3,lwd=.1)

    plot(function(t)ptnc(t, df = 10, ncp=35),0,70, log="y")
    ## "bug": jumps to 0 ..
    plot(function(t) pt(t, df = 10, ncp=35),0,70, add=TRUE, col="purple")
    abline(h=10^(-2:0),col="green",lty=3,lwd=.1)
}

p.pt <- function(args=alist(x, df=, ncp=), n=1001, log="")
{
  ## Purpose:
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  4 May 1999, 11:28

    p1 <-

    plot(function(x)ptnc(x, df=x, ncp=x), 1e-3, 1e3, n=1001, log="x")
    plot(function(x)  pt(x, df=x, ncp=x), 1e-3, 1e3, n=1001, add=TRUE,col=2)

}

par(mfrow=c(3,2),lab=c(10,20,7))
{
    plot(function(x)ptnc(x, df=x, ncp=x), 1e-3, 1e3, n=1001, log="x")
    plot(function(x)  pt(x, df=x, ncp=x), 1e-3, 1e3, n=1001, add=TRUE,col=2)
    ## two different jumps !
    ##
    plot(function(x)ptnc(x, df=x, ncp=2*x), 1e-3, 1e3, n=1001, log="x")
    plot(function(x)  pt(x, df=x, ncp=2*x), 1e-3, 1e3, n=1001, add=TRUE,col=2)
    ## look identical
    plot(function(x)ptnc(x, df=2*x, ncp=2*x), 1e-3, 1e3, n=1001, log="x")
    plot(function(x)  pt(x, df=2*x, ncp=2*x), 1e-3, 1e3, n=1001, add=TRUE,col=2)
}

###================= Central t  (ncp=0) ==== almost no differnce
## They overlap:
plot(function(t)ptnc(t,df=4,ncp=0),col="blue",-3,10)
plot(function(t)pt (t,df=4),add=T, col="red" ,-3,10)

## The small absolute Differences -- all negative ! --
plot(function(t)ptnc(t,df=4,ncp=0) - pt (t,df=4), col="blue",-3,10)
plot(function(t)ptnc(t,df=4,ncp=0) - pt (t,df=4), col="blue",0,1000,n=10001)
## The small relative Differences:
plot(function(t)1-pt (t,df=4)/ptnc(t,df=4,ncp=0),col="blue",-3,10)
plot(function(t)1-pt (t,df=4)/ptnc(t,df=4,ncp=0),col="blue",0,1000,n=10001)


## Compare dcdflib's ptnc() == Posten(1994)'s algorithm,
## with our pt() and with the  "simple approximation"s now from package 'DPQ'

## setwd("/u/maechler/R/MM/NUMERICS/dpq-functions")
## source("t-nonc-fn.R")## see also ./t-nonc-approx.R

## if(!dev.interactive(orNone=TRUE))
##    pdf("pt_ncp_appr.pdf", paper="default", height=0, width=0)

## FIXME: These do not even call ptnc() here !!!! <<<<<<<<<<<<<<<<<<
if(requireNamespace("DPQ")) {
    pntR    <- DPQ::pntR
    pntJW39 <- DPQ::pntJW39 ## <--
    op <- mult.fig(12, oma = c(0,0,3,1),
                   main = "pt(t, df=10, ncp) and approximations")$old.par
    ii <- seq(0,1, len=501)
    for(ncp in c(0,1,2,5,10,37.6218,37.6219,100,1000,1e4,1e5,1e10)) {
        tt <- ii * max(2,2.4*(ncp + if(ncp < 3) ncp else 0))
        print(c(system.time(p1 <- pt     (tt,df=10,ncp=ncp))[1],
                system.time(p2 <- pntR   (tt,df=10,ncp=ncp))[1],
                system.time(p.a<- pntJW39(tt,df=10,ncp=ncp))[1]))
        matplot(tt, cbind(p1,p2, p.a), col=c("black","blue3","red3"),
                lwd = c(1.5,1, .5), lty=1, type='l', xlab= "t", ylab = "",
                main = paste("ncp = ",formatC(ncp, digits=10,wid=1)))
        if(FALSE)
            abline(h=.1*(0:10),col='gray80', lwd=.3, lty=3)
        par(new=TRUE)
        matplot(tt, abs(cbind(p.a,p2) - p1), type="l",
                lty = 2:3, col = c("red3","blue3"), yaxt="n", ann = FALSE)
        if(ncp == 0)
            legend("topleft", c("|p.appr() - pt()|", "|ptnc() - pt()|"),
                   lty = 2:3, col = c("red3","blue3"),
                   inset=.02)
        axis(4, col="purple")
    }
    par(op)
}

## dev.off()

### End{ originally from .............../pnt-ex.R } -------------------------

p1 <- seq(0, 1, length= 11)
p2 <- seq(0, 1, length= 21)
p3 <- seq(0, 1, length=101)

unlist(pb <- cdfbet(which=1, p=0,q=0, x=.2,y=.8, a=2,b=3,0,0))
stopifnot(all.equal(pbeta(.2, 2,3), pb $ p, tol=1e-12))

## The new nice vectorized ones :
str(vp <- cdfbetV(which=1, p=0,q=0, x=.2,y=.8, a=2, b=2:20,0,0))
vp$p
stopifnot(all.equal(vp$p, pbeta(.2, 2, 2:20), tol=1e-12))

qbeta(p2, 2,3)

stopifnot(cdfbetV(which=2, p= p2, q= 1-p2, x=0,y=0, a=2,b=3,0,0)$x
          == qbetaB(p=p2, a=2, b=3)) # exactly : qbetaB() calls cdfbetV() !
all.equal(qbeta (p2, 2,3),
          qbetaB(p2, 2,3), tol=1e-10)
##"Mean relative difference: 3.831832e-09"


### Extreme tails:
cdfpoi(which=1, p=0,q=0, s=190, xlam = 100, status=0,b=0)$q
## 4.174239e-16
1-ppois(190,lam=100)
## 4.4408e-16

cdfpoi(which=1, p=0,q=0, s=195, xlam = 100, status=0,b=0)$q
## 1.4795033e-17
1-ppois(195,lam=100)
## 0 (!)

## cdfpoiV  --- now okay, once I use the RECYCLING *V functions:
ss <- seq(100, 195, by=5)
str(poi.ss <- cdfpoiV(which=1, p=0,q=0, s=ss, xlam = 100, status=0,b=0))
stopifnot(all.equal(poi.ss$p, ppois(ss, lam=100), tol = 1e-12),
          all.equal(poi.ss$q, 1 - poi.ss$p,       tol = 4e-15)
          )

##--- ptnc is ok
tt <- c(seq(-2,2,by=.25), 3:6)
p1 <- ptnc( tt, df=45, ncp=-.5)
p2 <- pt  ( tt, df=45, ncp=-.5)
stopifnot(all.equal(p1,p2, tol=1e-10))
summary(1 - p2/p1)# -4.295e-10 .. 4.794e-11

##--- This shows that   qtnc() now works :
qtnc(p= .975, df = 1:10, ncp=0)
##-  [1] 12.706205  4.302653  3.182446  2.776445  2.570582  2.446912  2.364624
##-  [8]  2.306004  2.262157  2.228139
str(tncV(which=2, p= .975, df = 1:10, ncp=0))

pp <- seq(.1,.9, by=.1)
stopifnot(all.equal(qtnc( p=pp , df=11, ncp=0),
                    qt  ( p=pp , df=11), tol = 2e-10)) # 1.09343e-10

qtnc( p=pp , df=11, ncp=1)
if (FALSE) ## Not allowed:
qt  ( p=pp , df=11, ncp=1)

## pbetaB()  from  ...../pq-funs.R :
x <- seq(0,1,len=101)
stopifnot(all.equal(pbetaB(x, .1,.2),
                    pbeta (x, .1,.2), tol = 1e-14),
          all.equal( pbetaB(x, .1, .2, lower.tail = FALSE),
                    1- pbeta(x, .1,.2), tol = 1e-14))

## Now (since May 18, when I changed `scale' for cdfgam), these match:
## --- to have back-compatibility: use  scale = 1 !!

### dcdflib.c (May 18) --- broken  gratio() !!!
##-- However, the following is ALL RIGHT for Original dcdflib-c !!!
xx <- seq(0,2, len=101); all.equal(pgamma (xx, 2), pgammaB(xx, 2),tol=0)
         ##Solaris> 2.667347 e-16   ## Intel Linux (2.0) 2.104258e-16
xx <- seq(2,15, len=101); all.equal(pgamma (xx, 2), pgammaB(xx, 2),tol=0)
         ##rel.diff:  Solaris> 1.2548 e-16   ## Intel Linux (2.0) 1.5896 e-16


## gratio() -- richtig : -- falsch im neuen dcdflib
print(gratio(2, .1,0,0, ind=0)$ans, dig=20)
##[1] 0.00467884016044447..

## grat1() -- dasselbe im neuen und alten dcdflib:
print(grat1(2, .1,0,0,0, eps=1e-20) $p, dig=20)
##[1] 0.00467884017669852..
print(grat1(2, .1,0,0,0, eps=1e-3) $p, dig=20)
##[1] 0.00467884017498421..

## Not good -- BAD for small x:
xx <- seq(0,2, len=101); all.equal(pgamma (xx, 2), pgammaB(xx, 2))# .00120!!
## Ok larger x
xx <- seq(2,15, len=101); all.equal(pgamma (xx, 2), pgammaB(xx, 2)) # TRUE

## This shows how bad it is:
xx <- seq(0,3, len=501)
plot(xx, pgamma (xx, 2) - pgammaB(xx, 2),
     main = "pgamma (xx, 2) - pgammaB(xx, 2)")
mtext(R.version.string, side=4, adj=0, cex = .75)
## print this ..

plot(xx, pgamma (xx, 1.5) - pgammaB(xx, 1.5))# still not perfect for x < a
plot(xx, pgamma (xx, .5) - pgammaB(xx, .5))# ok
## The next three are still not ok in new dcdflib.c ...
## Wrong only when compiled with -fpic (Linux 2.0; gcc 2.7.2.3)
plot(xx, pgamma (xx, .1) - pgammaB(xx, .1))# ok
plot(xx, pgamma (xx, 1e-4) - pgammaB(xx, 1e-4))# ok
plot(xx, pgamma (xx, 1e-14) - pgammaB(xx, 1e-14))# ok

xB <- 4^(1:80)
plot(xB, pgamma (xB, 1e-4) - pgammaB(xB, 1e-4), log='x')# all 0
xB <- 2^(1:80)
plot(xB, pgamma (xB, 1) - pgammaB(xB, 1), log='x')# all 0

## Not anymore all 0 --- "Asymptotic case":
##		(last 5 give exact 1/0, because r=0 after S40)
xB <- 1.05^(1:140)
plot(xB, pgamma (xB, 1.4) - pgammaB(xB, 1.4), log='x')
plot(xB, pgammaB(xB, 1.4), log='x')
plot(xB, pgammaB(xB, 1.4, lower.tail=FALSE), log='xy')

xB <- seq(100,200, len=201)
plot(xB, pgamma (xB, 144) - pgammaB(xB, 144), log='x')
plot(xB, pgammaB(xB, 144), log='x')
plot(xB, pgammaB(xB, 144, lower.tail=FALSE), log='xy')

## Gen. Temme Expansion: -- difference to pgamma() !
xB <- seq(1.2e4, 1.6e4, len=201)
plot(xB, pgamma (xB, 1.4e4) - pgammaB(xB, 1.4e4), log='x')
plot(xB, pgammaB(xB, 1.4e4), log='x')
plot  (xB, pgammaB(xB, 1.4e4, lower.tail= TRUE), log='xy')
points(xB, pgammaB(xB, 1.4e4, lower.tail=FALSE), col='blue')

## Gen. Temme Expansion: -- difference to pgamma() !
xB <- seq(1.398e7, 1.402e7, len=201)
plot(xB, pgamma (xB, 1.4e7) - pgammaB(xB, 1.4e7), log='x')
plot(xB, pgammaB(xB, 1.4e7), log='x')
plot  (xB, pgammaB(xB, 1.4e7, lower.tail= TRUE), log='xy')
points(xB, pgammaB(xB, 1.4e7, lower.tail=FALSE), col='blue')

## Gen. Temme Expansion: -- difference to pgamma() !
a <- 1.4e25
xB <- seq((1-1e-12)*a, (1+1e-12)*a, len=201)
plot(xB, pgamma (xB, a) - pgammaB(xB, a), log='x')
plot(xB, pgammaB(xB, a), log='x')
plot  (xB, pgammaB(xB, a, lower.tail= TRUE), log='xy')
points(xB, pgammaB(xB, a, lower.tail=FALSE), col='blue')



plot(function(x)pgamma (x, 2.5, 2.2), from=-1, to= 15)
plot(function(x)pgammaB(x, 1.5, 2.2), from=-1, to= 15, add=T, col=2)

plot(function(x)pgamma (x, 1.5, 2.2), from=-1, to= 15)
plot(function(x)pgammaB(x, 1.5, 2.2), from=-1, to= 15, add=T, col=2)


### qgammaB() <===> gaminv() testen !!!
qgammaB(p1,2,2)


### NOTA BENE: This has an "inverse" `scale' for cdfgamma() !
###			   ---------------      ==========

##==> at least  qtnc() work ok ! :
## qtnc() -> tncV() -> .C("V_cdftnc"..)
## ----      ----      ----------------
qtnc(p= .975, df = c(1:10,20,50,100), ncp=1)
##-  [1] 34.556218  8.626891  5.767894  4.792655  4.313083  4.030285  3.844446
##-  [8]  3.713235  3.615747  3.540506  3.229880  3.063452  3.010993
qtnc(p= .975, df = c(1:10,20,50,100), ncp=0)
##-  [1] 12.706205  4.302653  3.182446  2.776445  2.570582  2.446912  2.364624
##-  [8]  2.306004  2.262157  2.228139  2.085963  2.008559  1.983972

## The same from
str(cdftncV(which=2, p=.975, q=1-.975,t=pi,
            df=c(1:10,20,50,100), pnonc=0,
            status=pi,bound=0)[-(1:3)], vec.len=13)
