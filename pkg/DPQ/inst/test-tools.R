## Be *more* modular ==> do *not* load  Matrix test-tools here !
## source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))

## Just a small subset of those from 'Matrix' i.e,
##    system.file("test-tools-1.R",  package = "Matrix")

identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)
identical5 <- function(a,b,c,d,e) identical(a,b) && identical4(b,c,d,e)

assert.EQ <- function(target, current, tol = if(showOnly) 0 else 1e-15,
                      giveRE = FALSE, showOnly = FALSE, ...) {
    ## Purpose: check equality *and* show non-equality
    ## ----------------------------------------------------------------------
    ## showOnly: if TRUE, return (and hence typically print) all.equal(...)
    T <- isTRUE(ae <- all.equal(target, current, tolerance = tol, ...))
    if(showOnly) return(ae) else if(giveRE && T) { ## don't show if stop() later:
	ae0 <- if(tol == 0) ae else all.equal(target, current, tolerance = 0, ...)
	if(!isTRUE(ae0)) writeLines(ae0)
    }
    if(!T) stop("all.equal() |-> ", paste(ae, collapse=sprintf("%-19s","\n")))
    else if(giveRE) invisible(ae0)
}

pkgRversion <- function(pkgname)
    sub("^R ([0-9.]+).*", "\\1", packageDescription(pkgname)[["Built"]])

showSys.time <- function(expr, ...) {
    ## prepend 'Time' for R CMD Rdiff
    st <- system.time(expr, ...)
    writeLines(paste("Time", capture.output(print(st))))
    invisible(st)
}
showProc.time <- local({ ## function + 'pct' variable
    pct <- summary(proc.time())# length 3, shorter names
    function(final="\n", ind=TRUE) { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- summary(proc.time())
	delta <- (pct - ot)[ind]
	##  'Time' *not* to be translated:  tools::Rdiff() skips its lines!
	cat('Time', paste0("(",paste(names(delta),collapse=" "),"):"), delta, final)
    }
})

## We import from 'sfsmisc' anyway:
relErr  <- sfsmisc::relErr
relErrV <- sfsmisc::relErrV

##===================== DPQ-own ==================================================

## to be used in saveRDS(list_(nam1, nam2, ...),  file=*) :
list_ <- function(...) {
    ## nms <- vapply(sys.call()[-1L], deparse, "", width.cutoff=500L, backtick=FALSE)
    nms <- vapply(sys.call()[-1L], as.character, "")
    `names<-`(list(...), nms)
}
## even faster
list_ <- function(...)
   `names<-`(list(...), vapply(sys.call()[-1L], as.character, ""))

save2RDS <- function(x, file, do.time=TRUE, verbose=TRUE, ...) {
    if(verbose) cat("Saving to ", file, "\n")
    saveRDS(x, file=file, ...) # returning NULL
    if(do.time) showProc.time()
}
readRDS_ <- function(file, do.time=TRUE, verbose=TRUE, ...) {
    if(verbose) cat("Reading from ", file, "\n")
    if(do.time) on.exit(showProc.time())
    readRDS(file=file, ...)
}

##' load a named list `L` into environment `envir`
loadList <- function(L, envir = .GlobalEnv)
    invisible(lapply(names(L), function(nm) assign(nm, L[[nm]], envir=envir)))

## matplot(*, type = "b", ...)  but "smartly": not too many pch
matplotB <- function(x, y, type = "b", npts = 50L, ...) {
    stopifnot(type %in% c("p","o","b"), # the ones considered in matplot();  + "c" ?
              !is.null(dim(y)), length(npts) == 1L, is.numeric(npts), npts == as.integer(npts), npts >= 2L)
    x1 <- is.null(dx <- dim(x))
    n <- max(if(x1) length(x) else dx[1L], nrow(y)) # simple; other cases not covered
    if(n > npts) {
        matplot  (x, y, type="l", ...)
        i <- seq.int(1, n, length.out = npts)
        matpoints(if(x1) x[i] else x[i,], y[i,], type="p", ...)
    } else { # just
        matplot(x, y, type=type, ...)
    }
}

