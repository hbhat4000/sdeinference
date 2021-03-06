// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "usefnptr_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// callViaXPtr
arma::vec callViaXPtr(const arma::vec x, SEXP xpsexp);
RcppExport SEXP usefnptr_callViaXPtr(SEXP xSEXP, SEXP xpsexpSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xpsexp(xpsexpSEXP);
    __result = Rcpp::wrap(callViaXPtr(x, xpsexp));
    return __result;
END_RCPP
}
