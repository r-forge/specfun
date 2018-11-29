#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* Now .C() calls in ../R/toms644.R  --- C code in ./zbsubs.c */
#include "zbsubs.h"

static const R_CMethodDef CEntries[] = {
    {"zairy", (DL_FUNC) zairy,  8},
    {"zbesh", (DL_FUNC) zbesh, 10},
    {"zbesi", (DL_FUNC) zbesi,  9},
    {"zbesj", (DL_FUNC) zbesj,  9},
    {"zbesk", (DL_FUNC) zbesk,  9},
    {"zbesy", (DL_FUNC) zbesy, 11},
    {"zbiry", (DL_FUNC) zbiry,  8},
    {NULL, NULL, 0}
};

void R_init_Bessel(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
