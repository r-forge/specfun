#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran() calls in ../R/toms644.R  --- C code in ./zbsubs.c */
#include "zbsubs.h"

static const R_FortranMethodDef FortranEntries[] = {
    {"zairy", (DL_FUNC) &F77_NAME(zairy),  8},
    {"zbesh", (DL_FUNC) &F77_NAME(zbesh), 10},
    {"zbesi", (DL_FUNC) &F77_NAME(zbesi),  9},
    {"zbesj", (DL_FUNC) &F77_NAME(zbesj),  9},
    {"zbesk", (DL_FUNC) &F77_NAME(zbesk),  9},
    {"zbesy", (DL_FUNC) &F77_NAME(zbesy), 11},
    {"zbiry", (DL_FUNC) &F77_NAME(zbiry),  8},
    {NULL, NULL, 0}
};

void R_init_Bessel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
