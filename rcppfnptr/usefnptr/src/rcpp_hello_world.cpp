// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

#include "usefnptr_types.h"

// [[Rcpp::export]]
arma::vec callViaXPtr(const arma::vec x, SEXP xpsexp)
{
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  arma::vec y = fun(x);
  return(y);
}
 
