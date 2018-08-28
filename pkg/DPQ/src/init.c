#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

#include "DPQpkg.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


/* .C calls */

// qchisq_appr.c :
static R_NativePrimitiveArgType qchisq_appr_v_t[7] = {
    REALSXP, INTSXP, REALSXP, REALSXP,
    LGLSXP,  LGLSXP, REALSXP };
// pnchisq-it.c :
static R_NativePrimitiveArgType Pnchisq_it_t[] = {
    REALSXP, REALSXP, REALSXP,
    /* errmax: */ REALSXP, REALSXP, INTSXP,
    /* i_0, n_terms: */ INTSXP,  INTSXP,
    REALSXP, REALSXP };
// ppois-direct.c :
static R_NativePrimitiveArgType ppois_D_t[] = { REALSXP, INTSXP, REALSXP, REALSXP };

// wienergerm_nchisq.c :
static R_NativePrimitiveArgType pchisqV_t[] = {
    REALSXP, INTSXP,  REALSXP, REALSXP,
    LGLSXP,  LGLSXP,  INTSXP };



static const R_CMethodDef CEntries[] = {
    CDEF(qchisq_appr_v),
    CDEF(Pnchisq_it),
    CDEF(ppois_D),
    CDEF(pchisqV),

    {NULL, NULL, 0}
};

static R_FortranMethodDef FortEntries[] = {
    {"noncechi", (DL_FUNC) &F77_SUB(noncechi), 6}, // <- ./wienergerm_nchisq_F.f
    {NULL, NULL, 0}
};


void R_init_DPQ(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
