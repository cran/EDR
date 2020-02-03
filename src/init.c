#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>

void F77_NAME(edrstp3a)( int* d, int* n, int* dp1, double* wij, double* kksi,
  double* kksii, double* mat, double* wi, double* s, double* work, int* iwork);

void F77_NAME(edrstp3b)( int* d, int* n, int* dp1, double* wij, double* kksi,
  double* y, double* kksii, double* mat, double* s, double* u, double* vt,
  double* work, int* iwork, double* fx, double* fw, double* lll, double* yw);

void F77_NAME(edrstp3c)( int* d, int* n, int*nest, int* dp1, double* x,
  double* xest, double* kksi, double* y, double* mat, double* s, double* u,
  double* vt, double* work, int* iwork, double* fw, double* yw);

static R_NativePrimitiveArgType edrstp3a_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType edrstp3b_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType edrstp3c_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP};

static const R_FortranMethodDef fmethods[] = {
    						{"edrstp3a", (DL_FUNC) &edrstp3a_ , 11, edrstp3a_t},
                {"edrstp3b", (DL_FUNC) &edrstp3b_ , 17, edrstp3b_t},
                {"edrstp3c", (DL_FUNC) &edrstp3c_ , 16, edrstp3c_t},
                {NULL, NULL, 0, NULL}
};

void R_init_EDR(DllInfo *dll)
                {
                    R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
                    R_useDynamicSymbols(dll,FALSE);
                    R_forceSymbols(dll,TRUE);
                }
