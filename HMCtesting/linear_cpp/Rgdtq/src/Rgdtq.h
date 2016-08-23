#ifndef _RGDTQ_
#define _RGDTQ_

#include <RcppArmadillo.h>

RcppExport SEXP gdtqCPP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP dtqCPP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
RcppExport SEXP gdtqCPP_linear(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP dtqCPP_linear(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 

#endif

