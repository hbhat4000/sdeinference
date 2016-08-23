// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

arma::vec fun1_cpp(const arma::vec& x)
{
  arma::vec y = x + x;
  return (y);
}

arma::vec fun2_cpp(const arma::vec& x)
{
  arma::vec y = 10*x;
  return (y);
}

#include "Rdtq_types.h"

// [[Rcpp::export]]
XPtr<funcPtr> putFunPtrInXPtr(std::string fstr)
{
    if (fstr == "fun1")
        return(XPtr<funcPtr>(new funcPtr(&fun1_cpp)));
    else if (fstr == "fun2")
        return(XPtr<funcPtr>(new funcPtr(&fun2_cpp)));
    else
        return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

// [[Rcpp::export]]
arma::vec callViaString(const arma::vec x, std::string funname) {
    XPtr<funcPtr> xpfun = putFunPtrInXPtr(funname);
    funcPtr fun = *xpfun;
    arma::vec y = fun(x);
    return (y);
}

