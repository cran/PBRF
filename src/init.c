#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(ypbf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ypbff)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"ypbf",  (DL_FUNC) &F77_NAME(ypbf),  10},
  {"ypbff", (DL_FUNC) &F77_NAME(ypbff), 13},
  {NULL, NULL, 0}
};

void R_init_PBRF(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
