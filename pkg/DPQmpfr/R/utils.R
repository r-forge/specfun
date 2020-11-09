## Not exported, used to make  'R CMD check <pkg>'  be faster *or* more extensive:
doExtras <- function(int = interactive()) {
    int || nzchar(Sys.getenv("R_DPQmpfr_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

