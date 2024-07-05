### really related to
## vvvvvvvvvvvvvvvvvvvvvvvvvv
## ../../tests/stirlerr-tst.R   *and*
## ../../Misc/stirlerr/stirlerr-eval-tst.R
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  both source(<this>), i.e.
## source(system.file(package="DPQ", "extraR", "stirlerr-plots.R", mustWork=TRUE))

require("sfsmisc") # eaxis()

##' my R platform/architecture -- to be usable as file name (part) etc
myPlatform <- function(Rmin = 9L, osM = 12L)
    paste(abbreviate(sfsmisc::shortRversion(date=FALSE, spaces=FALSE), Rmin), # "R-devel..." too long
          .Platform$OS.type, # 'unix'
          sub("_$","", gsub("[^[:alnum:]]", "_",
                            abbreviate(osVersion, minlength=osM))),
          if(!capabilities("long.double") &&
             !grepl("aarch64-apple", R.version$platform)) "noLD",
          sep='_')

## "FIXME" other versions of mtextVersion() --> ~/R/Pkgs/Rmpfr/tests/special-fun-dgamma.R
if(FALSE) { #					~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  myRversion <- paste(R.version.string, "--", osVersion)
  if((mach <- Sys.info()[["machine"]]) != "x86_64")
      myRversion <- paste0(myRversion, "_", mach)
  if(!capabilities("long.double"))
      myRversion <- paste0(myRversion, "_no_LDbl")
  myRversion

  ## and then used as
  mtext(myRversion, adj=1, cex=3/4)
}
mtextVersion <- function(adj = 1, col = 1) {
    mtext(osVersion, line=1, col=col, adj=adj)
    mtext(sfsmisc::shortRversion(spaces=FALSE), col=col, adj=adj)
}

## Currently only used in >>>>>>>>>>>> ../../tests/stirlerr-tst.R  <<<<<<<<<<<<<<<<<<<<<<<
##                                     ~~~~~~~~~~~~~~~~~~~~~~~~~~
##                    lapply(setNames(,k.), function(k) find1cuts(k=k, c1=c1..$c1))
plot1cuts <- function(res1c, # resulting from
                      col = c(2,4),
                      ymin = max(min(y), 2e-17), ...)
{
    rnms <- c("k", "n", "relE", "smooths", "i.n")
    stopifnot(is.list(res1c), rnms %in% names(res1c))
    list2env(res1c, envir = environment())
    ## ---> { k, n, relE, smooths, i.n , ..} now exist in local environment
    stopifnot(exprs = {
        is.list(smooths)
        (nS <- length(smooths)) >= 1L
        sapply(smooths, dim) == rep(c(length(n), 2L), nS)
        identical(dim(i.n), c(2L, nS))
    })
    y <- abs19(relE)
    matplot(n, y, type = "l", log = "y", col=col, lwd=1/2, yaxt = "n", ylim = c(ymin, max(y, 5e-16)),
            xlab = quote(n), ylab = quote(abs(relE(n))),
            main = substitute({k == K} *";" ~~ n %in% group("[",list(N1,N2),"]"),
                              list(K = k, N1=signif(min(n),3), N2=signif(max(n),3))),
            ...)
    eaxis(2)
    drawEps.h()
    ## NB: all smoothing --- as in stirlerrPlot() above -- should happen in log-space
    ## (log(abs19(relEc)), 2, function(y) exp(smooth.spline(y, df=4)$y))
    ly <- log(y) # == log(abs19(relE)) == log(max(|r|, 1e-19))
    lines1 <- function(sy, ...) lines(n, sy, lwd=4, col=adjustcolor(col[1], 4/5), ...)
    lines2 <- function(sy, ...) lines(n, sy, lwd=4, col=adjustcolor(col[2], 1/2), ...)
    if(do.spl <- is.matrix(m <- smooths$spl)) {## add lines(smooth.spline())
        lines1(m[,"s1"])
        lines2(m[,"s2"])
    }
    if(do.low <- is.matrix(m <- smooths$low)) { ## lowess
        lines1(m[,"s1l"], lty=2)
        lines2(m[,"s2l"], lty=2)
    }
    ## also use cobs() splines for the 90% quantile !!
    if(do.cobs <- is.matrix(m <- smooths$cobs)) { ## lowess
        lines1(m[,"cs1"], lty=3)
        lines2(m[,"cs2"], lty=3)
    }
    had.n <- FALSE
    for(j in seq_len(ncol(i.n))) {
        ## i  <- i.n["i" ,j]  (unused)
        n. <- i.n["n.",j]
        if(length(n.) && is.finite(n.)) {
            if(!had.n) { # draw axis label only once
                had.n <- TRUE
                axis(3, at = signif(n.,3), col=col[1], line = -1)
            }
            abline(v = n., lty=2, lwd=3, col=adjustcolor(1, 1/2))
        }
    }
} ## plot1cuts()

## to enhance  |rel.Err| plots:  {also in ./pow-tst.R and  ~/R/Pkgs/Rmpfr/tests/special-fun-ex.R }
drawEps.h <- function(p2 = -(53:51), negative=FALSE, side = 4, lty = 3, lwd = 2, col = adjustcolor(2, 1/2)) {
    twop <- if(negative) c(outer(2^p2, c(-1,1))) else 2^p2
    labL <- lapply(p2, function(p) substitute(2^E, list(E=p)))
    if(negative) labL <- c(lapply(p2, function(p) substitute(-2^E, list(E=p))), labL)
    abline(h = twop, lty=lty, lwd=lwd, col=col)
    if(negative) abline(h=0, col=adjustcolor("gray20", 1/2), lwd=max(1, lwd))# at least thick as others
    axis(side, las=2, line=-1, at = twop, labels = as.expression(labL),
         col.axis = col, col=NA, col.ticks=NA)
}

