// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

arma::vec mycube(const arma::vec& x)
{
  arma::vec y = x % (x % x);
  return (y);
}

typedef arma::vec (*funcPtr)(const arma::vec& x);

// [[Rcpp::export]]
XPtr<funcPtr> mycubeXPtr()
{
  return(XPtr<funcPtr>(new funcPtr(&mycube)));
}

