## Be *more* modular ==> do *not* load  Matrix test-tools here !
## source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))

##--> showProc.time(), assertError(), relErrV(), ...
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

##' load a named list
loadList <- function(L, envir = .GlobalEnv)
    invisible(lapply(names(L), function(nm) assign(nm, L[[nm]], envir=envir)))
