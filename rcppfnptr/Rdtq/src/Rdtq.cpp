#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include "Rdtq_types.h"

// [[Rcpp::export]]
double callViaXPtr(const double x, SEXP xpsexp)
{
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  double y = fun(x);
  return(y);
}

// [[Rcpp::export]]
// some things that we need to add to this function:
// 1) ability to specify arbitrary domain boundaries and k
// 2) ability to specify xvec itself!
// 3) ability to return both xvec and the pdf
// 4) ability to make initial condition a pdf instead of a delta
NumericVector rdtq(double h, double k, int bigm, double init, double T, SEXP driftsexp, SEXP diffsexp)
{
  XPtr<funcPtr> driftF(driftsexp);
  XPtr<funcPtr> diffF(diffsexp);
  funcPtr driftfun = *driftF;
  funcPtr difffun = *diffF;

  unsigned int veclen = 2*bigm+1;

  NumericVector oldphatn(veclen);
  NumericVector phatn(veclen);
  NumericVector phatnp1(veclen);

  // pdf after one time step
  double mydiff = std::abs(difffun(init));
  double mydiff2 = pow(mydiff,2);
  for (int i=-bigm;i<=bigm;i++) {
    phatn(i+bigm) = exp(-pow(i*k-init-driftfun(init)*h,2)/(2.0*mydiff2*h))/(mydiff*sqrt(2.0*M_PI*h));
  }

  // iterate
  int bign = ceil(T/h);
  double thresh = 2.2e-16; // GSL_DBL_EPSILON;
  double lthresh = log(thresh);
  for (int n=1;n<bign;n++) {
    for (int i=-bigm;i<=bigm;i++) {
      bool keepgoing = true;
      int j = i;
      double tally = 0.0;
      while (keepgoing) {
        double thisdiff = pow(difffun(j*k),2);
        double lker = log(k) - pow(i*k - j*k - driftfun(j*k)*h,2)/(2.0*thisdiff*h) - 0.5*log(2.0*M_PI*thisdiff*h);
        if (lker < lthresh) keepgoing = false;
        if (phatn(j+bigm) >= thresh)
          tally += exp(lker + log(phatn(j+bigm)));
        j++;
        if (j > bigm) keepgoing = false;
      }
      if (i > -bigm) {
        keepgoing = true;
        j = i-1;
      }
      while (keepgoing) {
        double thisdiff = pow(difffun(j*k),2);
        double lker = log(k) - pow(i*k - j*k - driftfun(j*k)*h,2)/(2.0*thisdiff*h) - 0.5*log(2.0*M_PI*thisdiff*h);
        if (lker < lthresh) keepgoing = false;
        if (phatn(j+bigm) >= thresh)
          tally += exp(lker + log(phatn(j+bigm)));
        j--;
        if (j < -bigm) keepgoing = false;
      }
      phatnp1(i+bigm) = tally;
    }
    phatn = phatnp1;
  }
  return(phatn);
}



