#### Used to be part of /u/maechler/R/MM/NUMERICS/hyper-dist.R -- till Jan.24 2020

#### First version committed is the code "as in Apr 1999":
####  file | 11782 | Apr 22 1999 | hyper-dist.R


options(rErr.eps = 1e-30)
rErr <- function(approx, true, eps = .Options$rErr.eps)
{
    if(is.null(eps)) { eps <- 1e-30; options(rErr.eps = eps) }
    ifelse(abs(true) > eps,
           1 - approx / true, # rel.error
           true - approx) # absolute error in the limit
}

### ----------- -----------

k <- c(10*(1:9),100*(1:9),1000*(1:9))

options(digits=4)
ph.k2k <- phyper(k,2*k,2*k,2*k) # rel.err ~ 10^-7
cbind(k=k, phyper= ph.k2k,
      rERR.Ibeta= rErr(phyper.Ibeta     (k,2*k,2*k,2*k) , ph.k2k),
      rERR.as151= rErr(phyper.appr.as152(k,2*k,2*k,2*k) , ph.k2k),
      rERR.1mol = rErr(phyper.1molenaar (k,2*k,2*k,2*k) , ph.k2k),
      rERR.2mol = rErr(phyper.2molenaar (k,2*k,2*k,2*k) , ph.k2k),
      rERR.Peizer=rErr(phyper.Peizer    (k,2*k,2*k,2*k) , ph.k2k))

## Here, the Ibeta  fails (NaN) ; moleaar's are both very good :
ph.k2k <- phyper(k, 1.2*k, 2*k,1.5*k)
cbind(k=k, phyper= ph.k2k,
      rERR.Ibeta= rErr(phyper.Ibeta     (k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.as151= rErr(phyper.appr.as152(k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.1mol = rErr(phyper.1molenaar (k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.2mol = rErr(phyper.2molenaar (k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.Peizer=rErr(phyper.Peizer    (k,1.2*k,2*k,1.5*k) , ph.k2k))

x <- round(.8*k); ph.k2k <- phyper(x,1.6*k, 2*k,1.8*k)
(ph.mat <-
cbind(k=k, phyper= ph.k2k,
      rERR.Ibeta= rErr(phyper.Ibeta     (x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.as151= rErr(phyper.appr.as152(x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.1mol = rErr(phyper.1molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.2mol = rErr(phyper.2molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.Peiz = rErr(phyper.Peizer    (x,1.6*k,2*k,1.8*k) , ph.k2k))
)
matplot(k, abs(ph.mat[,-(1:2)]), type='o', log='xy')

par(mfrow=c(2,1))
for(Reps in c(0,1)) {
    options(rErr.eps = Reps) ## 0: always RELATIVE error; 1: always ABSOLUTE
    x <- round(.6*k); ph.k2k <- phyper(x,1.6*k, 2*k,1.8*k)
    print(ph.mat <-
          cbind(k=k, phyper= ph.k2k,
           rERR.Ibeta= rErr(phyper.Ibeta     (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.as151= rErr(phyper.appr.as152(x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.1mol = rErr(phyper.1molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.2mol = rErr(phyper.2molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.Peiz = rErr(phyper.Peizer    (x,1.6*k,2*k,1.8*k) , ph.k2k))
          )
    ## The two molenaar's  ``break down''; Peizer remains decent!
    Err.phyp <- abs(ph.mat[,-(1:2)])
    matplot(k, Err.phyp,
            ylim=
            if((ee <- .Options$rErr.eps) > 0) c(1e-200,1),
            type='o', log='xy')
}


hp1 <- list(n = 10, n1 = 12, n2 = 2)

for(h.nr in 1:1) {
  attach(hplist <- get(paste("hp",h.nr, sep='')) ) # defining  n, n1, n2
  cat("\n");  print(unlist(hplist))
  N <- n1 + n2
  x <- 0:max(0,min(n1, n2 - n))
  d <- dhyper(x, n1, n2, n)
  names(d) <- as.character(x)
  print(d)
}

### ----------- ----------- lfastchoose() etc -----------------------------


sapply(0:10, function(n) lfastchoose(n, c(0,n)))

## System  lchoose  gives non-sense:
sapply(0:10, function(n)    lchoose(n, c(0,n)))
## This is ok:
sapply(0:10, function(n) my.lchoose(n, c(0,n)))


###---------- some  gamma testing: -------------

pi - gamma(1/2)^2
for(n in 1:20) cat(n,":",formatC(log(prod(1:n)) - lgamma(n+1)),"\n")
for(n in 1:20) cat(n,":",formatC(prod(1:n) - gamma(n+1)),"\n")

## in  math/gamma.c :
p1 <- c(0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2)

### Bernoulli numbers  Bern() :

n <- 0:12
Bn <- n; for(i in n) Bn[i+1] <- Bern(i); cbind(n, Bn)
unix.time(Bern(30))
unix.time(Bern(30))
rm(.Bernoulli)
unix.time(Bern(30))
unix.time(Bern(80))

### as.lgamma() --- Asymptotic log gamma function :

as.lgamma(3,0)
for(n in 0:10) print(log(2) - as.lgamma(3,n))
for(n in 0:12) print((log(2*3*4*5) - as.lgamma(5+1,n)) / .Machine$double.eps)
for(n in 0:12) print((log(720)     - as.lgamma(6+1,n)) / .Machine$double.eps)
for(n in 0:12) print((log(720)     -    lgamma(6+1)  ) / .Machine$double.eps)
for(n in 0:12) print((log(5040)    - as.lgamma(7+1,n)) / .Machine$double.eps)


## look ok:
curve(lgamma(x),-70,10, n= 1001) # ok
curve(lgamma(x),-70,10, n=20001) #  slowness of x11 !!

curve(abs(gamma(x)),-70,10, n=  201, log = 'y')
curve(abs(gamma(x)),-70,70, n=10001, log = 'y')
par(new=T)
curve(lgamma(x),    -70,70, n= 1001, col = 'red')

## .., we should NOT  use log - scale with negative values... [graph ok, now]
curve(gamma(x),-     40,10, n=  2001, log = 'y')
curve(abs(gamma(x)),-40,10, n=  2001, log = 'y')
##- Warning: NAs produced in function "gamma"

x <- seq(-40,10, length = 201)
plot(x, gamma(x), log = 'y', type='h')


###------------------- early dhyper() problem ---------------------------------

##- Date: 30 Apr 97 11:05:00 EST
##- From: "Ennapadam Venkatraman" <VENKAT@biosta.mskcc.org>
##- Subject: R-beta: dhyper bug??
##- To: "r-testers" <r-testers@stat.math.ethz.ch>
##- Sender: owner-r-help@stat.math.ethz.ch

## I get the following incorrect answer  ---- FIXED

dhyper(34,410,312,49)
#   [1] 2.244973e-118
##- when the coorect value is
##-
dhyper(34,410,312,49)#          (from S-Plus)
##-    [1] 0.0218911

##================== R is now correct !
##E.S. Venkatraman (venkat@biosta.mskcc.org)


### phyperR() === R version of  {old} C version in  <R>/src/nmath/phyper.c :

## This takes long, currently {18 sec, on sophie [Ultra 1]}
k <- (1:100)*1000
unix.time(phyper.k <- phyper(k, 2*k, 2*k, 2*k))

k <- 1000; phyperR(k, 2*k, 2*k, 2*k) - phyper.k[1]

## Debug: ------------
if(interactive())
debug(phyperR)
phyperR(k, 2*k, 2*k, 2*k)
## before while():
xb <- 2000
xr <- 0
ltrm <- -2768.21584356709
NR <- 2000
NB <- 0# new value

