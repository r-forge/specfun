#include "DPQpkg.h" // before stdio (for MINGW_...)

#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


/* .C calls */

// qchisq_appr.c :
static R_NativePrimitiveArgType qchisq_appr_v_t[7] = {
    REALSXP, INTSXP, REALSXP, REALSXP,
    LGLSXP,  LGLSXP, REALSXP };
// pnchisq-it.c :
static R_NativePrimitiveArgType Pnchisq_it_t[] = {
    /* x, f, theta : */ REALSXP, REALSXP, REALSXP,
    /* errmax, reltol, itrmax, verbose: */ REALSXP, REALSXP, INTSXP, INTSXP,
    /* i_0, n_terms: */ INTSXP,  INTSXP,
    /* terms, prob : */ REALSXP, REALSXP };
static R_NativePrimitiveArgType ncbeta_t[] = {
    REALSXP, REALSXP, REALSXP, REALSXP,
    /* n: */ INTSXP,
    /* errmax: */ REALSXP, INTSXP, INTSXP,
    /* res: */ REALSXP };


// wienergerm_nchisq.c :
static R_NativePrimitiveArgType pchisqV_t[] = {
    REALSXP, INTSXP,  REALSXP, REALSXP,
    LGLSXP,  LGLSXP,  INTSXP };



static const R_CMethodDef CEntries[] = {
    CDEF(qchisq_appr_v),
    CDEF(Pnchisq_it),
    CDEF(ncbeta),
    CDEF(pchisqV),

    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
    CALLDEF(R_algdiv,   2), // <-- ./algdiv.c
    CALLDEF(ppoisD,     3), // <-- ./ppois-direct.c
    CALLDEF(Pnchisq_R, 14), // <-- ./pnchisq-it.c

    {NULL, NULL, 0}
};

static R_FortranMethodDef FortEntries[] = {
    {"noncechi", (DL_FUNC) &F77_SUB(noncechi), 6}, // <- ./wienergerm_nchisq_F.f
    {NULL, NULL, 0}
};


void R_init_DPQ(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
