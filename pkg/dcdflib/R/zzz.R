if(FALSE)
    ## library.dynam() not used when useDynLib(*) in ../NAMESPACE
.onLoad <- function(lib, pkg) {
    library.dynam("dcdflib", pkg, lib)
    ## TODO: use useDynLib(*) in NAMESPACE instead
}
