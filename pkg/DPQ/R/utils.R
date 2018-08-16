## Convenient argument checking

all_mpfr <- function(...) {
    for(x in list(...)) if(!inherits(x, "mpfr")) return(FALSE)
    ## else return
    TRUE
}

any_mpfr <- function(...) {
    for(x in list(...)) if(inherits(x, "mpfr")) return(TRUE)
    ## else return
    FALSE
}

