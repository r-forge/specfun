### TODO:  export  and document (ESS : C-c C-o C-r	ess-roxy-preview-Rd )
### ----   ------      --------

## NB: Currently explored (somewhat) in
## ---->  ~/R/Pkgs/DPQ/tests/bd0-tst.R
##          ~~~~~~~~~~~~~~~~~~~~~~~~~~

##' @title Vector of interesting x-values for dpois(x, lambda)
x4dpois <- function(lambda, f.max = 2) {
    if(lambda < 2000) 0:ceiling(max(f.max * lambda,
                                    qpois(1e-100, lambda, lower.tail=FALSE)))
    else if(lambda < 2^30) {
        s <- sqrt(lambda)
        seq(as.integer(lambda - (44-8e3/lambda)*s),
            as.integer(lambda + (44-8e3/lambda)*s), by=1L+ as.integer(lambda/25000))
    } else {# lambda >= 2^30 ~= 1.0737e9
        s <- max(lambda/1e14, sqrt(lambda))
        ## lambda/1e14: so that numerically lambda - 50*s < lambda !
        ## still want integer valued :
        seq(round(lambda - 50*s), round(lambda + 50*s), by=round(100*s / 1024))
    }
}

##' @title The "relative Error" in dpois()  --- using Rmpfr for ground truth
##' @param x, lambda, log  The arguments to dpoisFUN(), see \code{\link{dpois}}.  Note that \code{x}
##'                        gets a smart default derived from \code{lambda} and \code{f.maxX}.
##' @param prBits precision in bits to be used for "ground truth" \pkg{Rmpfr}'s \code{\link{dpois}()}.
##' @param f.maxX positive number, passed to \code{x4dpois}(*), the default for \code{x}.
##' @param errOnly  logical indicating if only `relE` the vector of relative errors should be returned.
##' @param keepMpfr logical indicating if the full Rmpfr() part should be returned as well.
##' @param dpoisFUN typically stats::dpois() but can use DPQ alternatives now, must have
##'         first three arguments equal \code{(x, lambda, log)}.
##' @param ...   optional further arguments to dpoisFUN(), notably if that's experimental
##'
##' @note  The result can be nicely graphed via   p.dpoisEr()
dpoisEr <- function(lambda, prBits = 1536, ## prBits = 256 is too small for lambda=1e100 !!!
                    log = FALSE
                  , f.maxX = 2
                  , x = x4dpois(lambda, f.max = f.maxX)
                  , errOnly = FALSE, keepMpfr = FALSE
                  , log2.min = .Machine$double.min.exp # .mpfr_erange("Emin") is *much* smaller!
                  , dpoisFUN = stats::dpois, ...)
{
    dpoisNam <- deparse(substitute(dpoisFUN))
    d... <- if(...length()) ## the  ` ... `  (as string)
                sub("^pairlist\\(","",
                    sub("\\)$","",
                        deparse1(match.call(expand.dots=FALSE)[["..."]])))
    stopifnot(is.numeric(prBits), length(prBits) == 1, prBits >= 64, as.integer(prBits) == prBits,
              is.logical(log), length(log) == 1, !is.na(log), length(lambda) == 1,
              is.function(dpoisFUN), length(fArgs <- formals(dpoisFUN)) >= 3,
              names(fArgs)[1:3] == c("x", "lambda", "log"))
    dpR <- dpoisFUN    (x, lambda, log=log, ...)
    dpM <- Rmpfr::dpois(x, lambda = mpfr(lambda, precBits=prBits), log=log)
    ## relative error, assuming dpM is  "the truth" :
    relE <- asNumeric(dpR/dpM - 1)
    ## must use *absolute* error when "truth" is zero:
    ## NOTE: mpfr will underflow *much* later than double.xmin:
    ##      mpfr(2,22)^.mpfr_erange("Emin")  # |--> 4.7651298e-323228497
    if(had0 <- any(d0 <- log2(abs(dpM)) < log2.min))
        ## dpR "must have underflown" here: replace rel.error values by absolute ones:
        relE[d0] <- asNumeric(dpR[d0] - dpM[d0])
                                        # the latter is typically 0 (as relE[.] was NaN)
    if(errOnly)
        relE
    else { ## return more {but via *attributes*, so the result can still be treated like a numeric}:
        structure(relE,
                  x=x, dpoisNam = dpoisNam, dpR=dpR, lambda=lambda, log=log,
                  had0 = had0, i0 = if(had0) which(d0), # where did we return *absolute* errors
                  call = match.call(),
                  d... = d...,
                  dpMpfr = if(keepMpfr) dpM)
    }
} ## {dpoisEr}


##' Deparse a *seq*uence, trying  seq() form (for arithmetic seq.)
deparse1seq <- function(x, nonSmall = 37,
                        ##            ^^ (large but not for plot here)
                        digits = 6,
                        ##       ^ using = 10 .. hmm default ?
                        wEnd = nonSmall %/% 4, epsD = 1e-5) {
    r <- deparse1(x) # in integer cases, already gives  <n>:<m>  {but not in "double"!}
    if((w <- nchar(r)) >= nonSmall) { # try arithmetic sequence  seq(a, z, by = delta)
	if(length(ud <- unique(diff(x))) == 1 ||
           max(abs(ud/(ud <- mean(ud)) - 1)) < epsD)
        {
            n <- length(x); xr <- x[c(1,n)] # xr = x-range
            f0 <- if(all(abs(xr) < .Machine$integer.max) &&
                     all(xr == as.integer(xr))) "%d"
                  else if(digits == 6) "%g"
                  else paste0("%.",digits,"g")
            f.by <- if(f0 == "%d") "%g" else f0
            fmt <- sprintf("seq(%s, %s, by=%s)", f0, f0, f.by)
	    sprintf(fmt, x[1], x[length(x)], ud)
        } else
	    paste0(substr(r, 1, nonSmall- wEnd - 4), "... ..", substr(r, w-wEnd, w))
    }
    else
	r
}

##' format() a (min,max)  range()-alike into  (closed) interval notation
formatRng <- function(rng, digits=4, newLine=TRUE, closed=c(TRUE,TRUE), ...) {
    stopifnot(is.logical(closed), !anyNA(closed), (lc <- length(closed)) >= 1)
    if(lc == 1) closed <- rep_len(closed, 2) # recycling
    paste0(if(closed[1]) "[" else "(",
           paste(format(rng, digits=digits, ...), collapse=", "),
           if(closed[2]) "]" else ")",
           if(newLine)"\n")
}

##' plot()  dpoisEr() explorations nicely
##' @param rEstr list, resulting from dpoisEr() calls
##' @param cex
##' @param xlab, main x-axis label and main title
##' @param col
##' @param lwd2, col2 : linewidth and color, needed \code{if(showP)}.
##' @param showP logical specifying if ...
##' @param epsLine logical indicating if dotted horizontal lines at +- eps[C] should be drawn.
##' @param col.eps color for eps[C]-lines
##' @param smallMid
##' @param epsSmall
##' @param small.bd0 positive number or FALSE, specifying the \emph{range} of small |x-lambda|
##'           values should be vizualized (by two vertical dashed lines in \code{col.bd0} color).
##' @param col.bd0
##' @param col.mid
##' @param ...
p.dpoisEr <- function(rEstr, cex = 1/4, xlab = deparse1seq(x, digits=10),
                     col = if(!had0) 1 else local({ cc <- rep(1, n); cc[at[["i0"]]] <- 2; cc }),
                     lwd2 = 4, col2 = adjustcolor(4, 1/4),
                     ## every showP-th (x, Pr(x)) value is shown; FALSE or 0: show nothing
                     showP = if(n <= 128) TRUE else round(n/100),
                     epsLine = TRUE, col.eps = adjustcolor("orange",1/2),
                     smallMid = TRUE,
                     epsSmall = if(log) (3+pmin(lambda,1e4)/1000)*1e-16
                                else    (2+pmin(lambda,1e4)/1000)*1e-14,
                     small.bd0 = if(isTRUE(grepl("ebd0", at[["call"]][["version"]]))) 0
                                 else if(length(del <- at[["call"]][["bd0.delta"]])) del else 0.1,
                     col.bd0 = 2, col.mid = 3,
                     main = sprintf("Relative Error of %s(x, lambda=%s%s%s)"
                                  , if(length(FN <- at[["dpoisNam"]])) as.character(FN) else "dpois"
                                  , lambda, if(log) ", log=TRUE" else ""
                                  , if(!is.null(.d <- at[["d..."]])) paste(",",.d) else ""
                                    ),
                     ...)
{
    stopifnot(is.numeric(rEstr), length(at <- attributes(rEstr)) >= 4,
              (n <- length(x <- at[["x"]])) == length(rEstr),
              is.numeric(lambda <- at[["lambda"]]), lambda > 0, length(lambda) == 1,
              is.logical(log    <- at[["log"]]),
              is.logical(had0   <- at[["had0"]]))
    ## a col *vector* is only applied to points, not lines ==> cannot use type = "b"
    plot (x, rEstr, cex=cex, col=col, main=main, xlab=xlab, ...)
    lines(x, rEstr,          col=adjustcolor(1, 1/3))
    abline(h= 0, lty=3, col = adjustcolor("gray40", 1/2))
    if(epsLine) {
        eps2 <- c(-1,1)*.Machine$double.eps
        abline(h= eps2, col=col.eps, lty=2)
        axis(2, at=eps2, labels = expression(-epsilon[C], +epsilon[C]),
             col=NA, col.lab=col.eps, col.axis=col.eps,
             las=2, line=-2.8, padj=c(1,0))
    }
    ## bd0() has had the "historical" cutoff  abs(x-np) < 0.1*(x+np)
    ## -- w/ np := lambda for dpois_raw()
    if(small.bd0 > 0) {
        d <- small.bd0 # = delta = 0.1 was hardwired in bd0()
        bd0Rng <- c((1-d)/(1+d), (1+d)/(1-d)) * lambda
        abline(v = bd0Rng, col=adjustcolor(col.bd0, 3/4), lty=2, lwd=3)
        u <- par("usr")
        text(max(u[1], bd0Rng[1]), u[4], paste("bd0.delta =", formatC(d)),
             col=col.bd0, adj = c(-.1, 1.1))
        cat("bd0_range for (d=",formatC(d),"): ", formatRng(bd0Rng), sep="")
    }
    if(smallMid) {
        sml.x <- abs(rEstr) < epsSmall
        sml.mid <- as.logical(c(
            rev(cumprod(rev(sml.x[x <  lambda]))),
            cumprod(        sml.x[x >= lambda])))
        if(any(sml.mid)) {
            midRng <- range(x[sml.mid])
            cat("middle (rE < epsS=",formatC(epsSmall),"): ", formatRng(midRng), sep="")
            points(rEstr ~ x, subset = sml.mid, col = col.mid, pch=".", cex=1/2)
            abline(v = midRng, col=adjustcolor(col.mid, 1/2))
            axis(1, at=midRng, col=NA, col.axis=col.mid, lwd=2, line=-2)
            mtext(substitute(epsilon[S] == EE, list(EE=formatC(epsSmall))),
                  side=3, at=lambda, col=col.mid, line=0)
        }
        else cat("empty middle (sml.mid is all FALSE)\n")
    }
    ## if(grepl("devel", R.version$status)) # on R-devel
    if(nzchar(R.version$status)) # on R-devel, R-patched ..
        mtext(R.version.string, 4, adj=0, cex=3/4)
    op <- par(new=TRUE); on.exit(par(op))
    if(showP) {
        if(is.null(dpM <- at[["dpMpfr"]])) {
            message("No 'dpMpfr' component, using 'dpR' instead; set keepMpfr=TRUE to change")
            dpM <- at[["dpR"]]
        }
        i <- if(showP == 1) TRUE else seq(1, n, by = showP)
        plot(x[i], asNumeric(dpM[i]), type="h",
             ann=FALSE, axes=FALSE, lwd=lwd2, col=col2)
    }
} ## p.dpoisEr()

##' Tabulate  dpoisEr() results
##' @param rEstr list, resulting from dpoisEr() calls
rE2mat <- function(rEstr) {
    stopifnot(length(dpM <- attr(rEstr,"dpMpfr")) > 0,
              is.numeric(dp <- as.vector(attr(rEstr, "dpR"))), length(dp) == length(dpM))
    cbind(x_m = drop(scale(attr(rEstr,"x"), center=TRUE,scale=FALSE)),
          dp = dp,
          dpM    = asNumeric(dpM),
          lg10dM = asNumeric(log10(dpM)),
          relE= asNumeric(relErrV(dp, target = dpM)))
}
