###
### Produce (reproducible) output --- for library(dcdflib) routines
###
library(dcdflib)

df1 <- c(1e-5,1e-2,.1, 2:10,20,50,100,1e3,1e5,1e20,1e100)

## tt could depend on ncp & df  (should really depend on  ncp !)
t0 <- c(-20,-1,-.2, 0, .2, .5, 1,2, pi, 10)
cat("\n\nptnc( ncp +t0, df, ncp), where t0 =", formatC(t0),"\n")
for(df in df1) {
    cat("\n---------\ndf=", formatC(df,wid=5),"\n=====\n")
    for(ncp in c(0, .1, 1, 10, 100, 1e10)) {
        tt <- t0 + ncp
        cat(" ncp=", formatC(ncp,wid=5),":",
            formatC(ptnc(tt, df=df, ncp=ncp),dig=7),"\n")
    }
}



## qtnc() looks really broken to me ! --- yes, it is also with v1.1 of dcdflib()
##
## Problem maybe in  tchiV ?

## MM (2018): Get "Warning: error flag ...." exactly for the first two df = 1e-5, 1e-2:
p0 <- c(1e-10,1e-3,.01, .1, .2, .5, .8, .9, .99)
cat("\n\nqtnc(p0, df, ncp), where p0 =", formatC(p0),"\n")
for(df in df1) {
    cat("\n---------\ndf=", formatC(df,wid=5),"\n=====\n")
    for(ncp in c(0, .1, 1, 10, 100, 1e10)) {
        cat(" ncp=", formatC(ncp,wid=5),":",
            formatC(qtnc(p0, df=df, ncp=ncp),dig=7),"\n")
    }
}


