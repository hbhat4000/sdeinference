#include <iostream>
#include <armadillo>
#include <cmath>
#include <gsl/gsl_math.h>

using namespace std;
using namespace arma;

double driftfun(double y)
{
    return(-y);
}

double difffun(double y)
{
    return(1.0);
}

int main(void) {
  double h = 0.001;
  double s = 0.75;
  double k = pow(h,s);
  double yM = datum::pi/k;
  int bigm = ceil(yM/k);
  unsigned int veclen = 2*bigm+1;
  double init = 0;

  vec oldphatn = zeros<vec>(veclen);
  vec phatn = zeros<vec>(veclen);
  vec phatnp1 = zeros<vec>(veclen);

  // pdf after one time step
  double mydiff = abs(difffun(init));
  double mydiff2 = pow(mydiff,2);
  for (int i=-bigm;i<=bigm;i++) {
    phatn(i+bigm) = exp(-pow(i*k-init-driftfun(init)*h,2)/(2.0*mydiff2*h))/(mydiff*sqrt(2.0*datum::pi*h));
  }

  // iterate
  double T = 1.0;
  int bign = ceil(T/h);
  double thresh = GSL_DBL_EPSILON;
  double lthresh = log(thresh);
  for (int n=1;n<bign;n++) {
    cout << "n = " << n << "\n";
    cout.flush();

    omp_set_dynamic(false);
    omp_set_num_threads(16);

#pragma omp parallel for

    for (int i=-bigm;i<=bigm;i++) {
      bool keepgoing = true;
      int j = i;
      double tally = 0.0;
      while (keepgoing) {
	double thisdiff = pow(difffun(j*k),2);
	double lker = log(k) - pow(i*k - j*k - driftfun(j*k)*h,2)/(2.0*thisdiff*h) - 0.5*log(2.0*datum::pi*thisdiff*h);
	if (lker < lthresh) keepgoing = false;
        if (phatn(j+bigm) >= thresh)
	  tally += exp(lker + log(phatn(j+bigm)));
	// cout << "i = " << i << ", j = " << j << ", ker = " << ker << '\n';
	j++;
	if (j > bigm) keepgoing = false;
      }
      if (i > -bigm) {
	keepgoing = true;
	j = i-1;
      }
      while (keepgoing) {
	double thisdiff = pow(difffun(j*k),2);
	double lker = log(k) - pow(i*k - j*k - driftfun(j*k)*h,2)/(2.0*thisdiff*h) - 0.5*log(2.0*datum::pi*thisdiff*h);
	if (lker < lthresh) keepgoing = false;
        if (phatn(j+bigm) >= thresh)
	  tally += exp(lker + log(phatn(j+bigm)));
	// cout << "i = " << i << ", j = " << j << ", ker = " << ker << '\n';
	j--;
	if (j < -bigm) keepgoing = false;
      }
      phatnp1(i+bigm) = tally;
    }

#pragma omp barrier

    phatn = phatnp1;
  }
  cout << "\n\n" << "********" << "\n\n" << phatn << "\n";
}  


