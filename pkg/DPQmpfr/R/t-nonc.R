## dntM() := version for Rmpfr:
nB <- length(body(dntM <- dnt))
body(dntM)[[nB]] <- substitute(new("mpfr", THIS), list(THIS = body(dnt)[[nB]]))
rm(nB) #                       -----------     -
