#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

#include "cdflib.h"
#include "vec_cdf.h"


#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


/* .C calls */

// a "list" of argument list types (that often appear)
#ifdef _not_yet_used_
static R_NativePrimitiveArgType C_i_t[1]  = { INTSXP };
static R_NativePrimitiveArgType C_ii_t[2] = { INTSXP, INTSXP };

static R_NativePrimitiveArgType C_d_t[1]  = { REALSXP };
static R_NativePrimitiveArgType C_dd_t[2] = { REALSXP, REALSXP };
#endif
static R_NativePrimitiveArgType C_d3_t[3] = { REALSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType C_d4_t[4] = { REALSXP, REALSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType C_d5_t[5] = { REALSXP, REALSXP, REALSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType C_d6_t[6] = { REALSXP, REALSXP, REALSXP, REALSXP,
					      REALSXP, REALSXP };

static R_NativePrimitiveArgType C_d4i_t[5] = { REALSXP, REALSXP, REALSXP, REALSXP, INTSXP };
static R_NativePrimitiveArgType C_d5i_t[6] = { REALSXP, REALSXP, REALSXP,
					       REALSXP, REALSXP, INTSXP };
static R_NativePrimitiveArgType C_d6i_t[7] = { REALSXP, REALSXP, REALSXP,
					       REALSXP, REALSXP, REALSXP, INTSXP };

#ifdef _not_yet_used_
static R_NativePrimitiveArgType C_id_t[2]  = { INTSXP, REALSXP };
static R_NativePrimitiveArgType C_idd_t[3] = { INTSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType C_id3_t[4] = { INTSXP, REALSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType C_id4_t[5] = { INTSXP, REALSXP, REALSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType C_id5_t[6] = { INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP };
#endif

static R_NativePrimitiveArgType C_id4id_t[7] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, REALSXP };
static R_NativePrimitiveArgType C_id5id_t[8] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, REALSXP };
static R_NativePrimitiveArgType C_id6id_t[9] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, REALSXP };


// in vec_cdf.h : -------------------------------------------------
static R_NativePrimitiveArgType C_id6idi_t[10] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, REALSXP, INTSXP };
static R_NativePrimitiveArgType C_id5idi_t[9] = { INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
						INTSXP, REALSXP, INTSXP };
static R_NativePrimitiveArgType C_id4idi_t[8] = { INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
						INTSXP, REALSXP, INTSXP };
#ifdef _not_yet_used_
static R_NativePrimitiveArgType C_id3idi_t[7] = { INTSXP, REALSXP, REALSXP, REALSXP,
						INTSXP, REALSXP, INTSXP };
#endif

static const R_CMethodDef CEntries[] = {
    {"bgrat",    (DL_FUNC) &bgrat,     7, C_d6i_t},
    {"bratio",   (DL_FUNC) &bratio,    7, C_d6i_t},

    {"cdfbet",   (DL_FUNC) &cdfbet,    9, C_id6id_t},
    {"cdfbin",   (DL_FUNC) &cdfbin,    9, C_id6id_t},
    {"cdfchi",   (DL_FUNC) &cdfchi,    7, C_id4id_t},
    {"cdfchn",   (DL_FUNC) &cdfchn,    8, C_id5id_t},
    {"cdff",     (DL_FUNC) &cdff,      8, C_id5id_t},
    {"cdffnc",   (DL_FUNC) &cdffnc,    9, C_id6id_t},
    {"cdfgam",   (DL_FUNC) &cdfgam,    8, C_id5id_t},
    {"cdfnbn",   (DL_FUNC) &cdfnbn,    9, C_id6id_t},
    {"cdfnor",   (DL_FUNC) &cdfnor,    8, C_id5id_t},
    {"cdfpoi",   (DL_FUNC) &cdfpoi,    7, C_id4id_t},
    {"cdft",     (DL_FUNC) &cdft,      7, C_id4id_t},
    {"cdftnc",   (DL_FUNC) &cdftnc,    8, C_id5id_t},

    {"cumbet",   (DL_FUNC) &cumbet,    6, C_d6_t},
    {"cumbin",   (DL_FUNC) &cumbin,    6, C_d6_t},
    {"cumchi",   (DL_FUNC) &cumchi,    4, C_d4_t},
    {"cumchn",   (DL_FUNC) &cumchn,    5, C_d5_t},
    {"cumf",     (DL_FUNC) &cumf,      5, C_d5_t},
    {"cumfnc",   (DL_FUNC) &cumfnc,    6, C_d6_t},
    {"cumgam",   (DL_FUNC) &cumgam,    4, C_d4_t},
    {"cumnbn",   (DL_FUNC) &cumnbn,    6, C_d6_t},
    {"cumnor",   (DL_FUNC) &cumnor,    3, C_d3_t},
    {"cumpoi",   (DL_FUNC) &cumpoi,    4, C_d4_t},
    {"cumt",     (DL_FUNC) &cumt,      4, C_d4_t},
    {"cumtnc",   (DL_FUNC) &cumtnc,    5, C_d5_t},

    {"gaminv",   (DL_FUNC) &gaminv,    6, C_d5i_t},
    {"grat1",    (DL_FUNC) &grat1,     6, C_d6_t},
    {"gratio",   (DL_FUNC) &gratio,    5, C_d4i_t},
    {"pni",      (DL_FUNC) &pni,       5, C_d4i_t},
/* vec_cdf.h : */
    {"V_cdfbet", (DL_FUNC) &V_cdfbet, 10, C_id6idi_t},
    {"V_cdfbin", (DL_FUNC) &V_cdfbin, 10, C_id6idi_t},
    {"V_cdfchi", (DL_FUNC) &V_cdfchi,  8, C_id4idi_t},
    {"V_cdfchn", (DL_FUNC) &V_cdfchn,  9, C_id5idi_t},
    {"V_cdff",   (DL_FUNC) &V_cdff,    9, C_id5idi_t},
    {"V_cdffnc", (DL_FUNC) &V_cdffnc, 10, C_id6idi_t},
    {"V_cdfgam", (DL_FUNC) &V_cdfgam,  9, C_id5idi_t},
    {"V_cdfnbn", (DL_FUNC) &V_cdfnbn, 10, C_id6idi_t},
    {"V_cdfnor", (DL_FUNC) &V_cdfnor,  9, C_id5idi_t},
    {"V_cdfpoi", (DL_FUNC) &V_cdfpoi,  8, C_id4idi_t},
    {"V_cdft",   (DL_FUNC) &V_cdft,    8, C_id4idi_t},
    {"V_cdftnc", (DL_FUNC) &V_cdftnc,  9, C_id5idi_t},

    {NULL, NULL, 0}
};

void R_init_dcdflib(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
