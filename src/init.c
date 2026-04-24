#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#define GFA_USING_R
#include "gfa.h"

/* Forward declarations of .Call() entry points defined in R_interface.c */
extern SEXP gfa_find_ir(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gfa_find_mr(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gfa_find_dr(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gfa_find_gq(SEXP, SEXP, SEXP);
extern SEXP gfa_find_zdna(SEXP, SEXP, SEXP);
extern SEXP gfa_find_str(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gfa_find_apr(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"gfa_find_ir",   (DL_FUNC)&gfa_find_ir,   7},
    {"gfa_find_mr",   (DL_FUNC)&gfa_find_mr,   5},
    {"gfa_find_dr",   (DL_FUNC)&gfa_find_dr,   5},
    {"gfa_find_gq",   (DL_FUNC)&gfa_find_gq,   3},
    {"gfa_find_zdna", (DL_FUNC)&gfa_find_zdna, 3},
    {"gfa_find_str",  (DL_FUNC)&gfa_find_str,  5},
    {"gfa_find_apr",  (DL_FUNC)&gfa_find_apr,  4},
    {NULL, NULL, 0}
};

void R_init_nonbgfa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
