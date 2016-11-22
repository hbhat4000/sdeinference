#include "Rdtq2d.h"

#include <RcppArmadillo.h>
#include <cmath>

using namespace std;
using namespace arma;

// gaussian pdf at the first step
static inline vec gaussian_pdf_vec_first(const vec& x, const double mu, const double sigma)
{
  vec u = (x - mu) / fabs(sigma);
  vec p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u % u) / 2);
  return p;
}

// point-wise gaussian pdfs for the internal steps
static inline double gaussian_pdf(const double x, const double mu, const double sigma)
{
  double u = (x - mu) / fabs(sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u * u) / 2);
  return p;
}

// yet another gaussian pdf, to use instead of interpolation at the final step
static inline vec gaussian_pdf_vec_last(const double x, const vec& mu, const vec& sigma)
{
  vec u = (x - mu) / abs(sigma);
  vec p = (1 / (sqrt (2 * M_PI) * abs (sigma))) % exp (-(u % u) / 2);
  return p;
}

// f function for 2d
static inline double f1(double x1, double x2, const vec& thetavec)
{
  return(-thetavec(0)*thetavec(0)*x2);
}

static inline double f2(double x1, double x2, const vec& thetavec)
{
  return(thetavec(1)*thetavec(1)*x1 + thetavec(2)*thetavec(2)*x1*x1*x1);
}

// g function for 2d
static inline double g1(double x1, double x2, const vec& thetavec)
{
  return(thetavec(3)*thetavec(3));
}

static inline double g2(double x1, double x2, const vec& thetavec)
{
  return(thetavec(4)*thetavec(4));
}

vec dtq(const vec &thetavec, const vec &C1, const vec &C2, double h, int numsteps, double k, double yM)
{
  double h12 = sqrt(h);
  int M = ceil(yM/k);
  int veclen = 2*M+1;
  int nr = veclen*veclen;
  int datapoints = C1.n_elem;
  vec xvec = k*linspace<vec>(-M,M,veclen);
  mat approxpdfvec(nr,datapoints-1);
  vec outpdf = zeros<vec>(datapoints-1);

// #ifdef _OPENMP
  omp_set_num_threads(24);
// #endif

if(numsteps > 1) { 
  #pragma omp parallel for
  for (int i=0; i<(datapoints-1); i++) {
    vec v1 = gaussian_pdf_vec_first(xvec, C1(i) + f1(C1(i),C2(i),thetavec)*h, h12*g1(C1(i),C2(i),thetavec));
    vec v2 = gaussian_pdf_vec_first(xvec, C2(i) + f2(C1(i),C2(i),thetavec)*h, h12*g2(C1(i),C2(i),thetavec));
    approxpdfvec.col(i) = vectorise(kron(v1, v2.t()));
  }
  #pragma omp barrier
}
else {
  #pragma omp parallel for
  for (int i=0; i<(datapoints-1); i++) {
    double v1 = gaussian_pdf(C1(i+1), C1(i) + f1(C1(i),C2(i),thetavec)*h, h12*g1(C1(i),C2(i),thetavec));
    double v2 = gaussian_pdf(C2(i+1), C2(i) + f2(C1(i),C2(i),thetavec)*h, h12*g2(C1(i),C2(i),thetavec));
    outpdf(i) = v1*v2;
  }
  #pragma omp barrier
  return outpdf;
}

  double supg = 0.5;
  int gamma = ceil(5*h12*supg/k);

  // loop over the timesteps
  for (int step=1; step<(numsteps-1); step++)
  {
    mat approxpdfvecnew = zeros<mat>(nr,datapoints-1);
#pragma omp parallel for
    for (int newspace=0; newspace<nr; newspace++)
    {
      int j = floor(newspace/veclen); // subtract M to get math i, j \in [-M,M]
      int i = newspace - veclen*floor(newspace/veclen);      
      for (int ip=(i-gamma); ip<=(i+gamma); ip++)
      {
        for (int jp=(j-gamma); jp<=(j+gamma); jp++)
        {
          int test = jp*veclen + ip;
          if ((test >= 0) && (test < nr))
          {
            int ipeff = ip - i + gamma;
            int jpeff = jp - j + gamma;
            if ((ipeff >= 0) && (ipeff <= 2*gamma) && (jpeff >= 0) && (jpeff <= 2*gamma))
            {
              double y1 = (ip-M)*k;
              double y2 = (jp-M)*k;
              double mu1 = y1 + f1(y1,y2,thetavec)*h;
              double mu2 = y2 + f2(y1,y2,thetavec)*h;
              double sigma1 = h12*g1(y1,y2,thetavec);
              double sigma2 = h12*g2(y1,y2,thetavec);
              double locbigg1 = gaussian_pdf(xvec(i), mu1, sigma1);
              double locbigg2 = gaussian_pdf(xvec(j), mu2, sigma2);
              approxpdfvecnew.row(newspace) += k*k*locbigg1*locbigg2*approxpdfvec.row(test);
            }
          }
        }
      }
    }
#pragma omp barrier
    approxpdfvec = approxpdfvecnew;
  }

  vec mu1 = zeros<vec>(nr);
  vec mu2 = zeros<vec>(nr);
  vec sig1 = zeros<vec>(nr);
  vec sig2 = zeros<vec>(nr);

#pragma omp parallel for
  for (int r=0; r<nr; r++)
  {
    int j = floor(r/veclen); // subtract M to get math i, j \in [-M,M]
    int i = r - veclen*floor(r/veclen);      
    mu1(r) = xvec(i) + f1(xvec(i),xvec(j),thetavec)*h;
    mu2(r) = xvec(j) + f2(xvec(i),xvec(j),thetavec)*h;
    sig1(r) = h12*g1(xvec(i),xvec(j),thetavec);
    sig2(r) = h12*g2(xvec(i),xvec(j),thetavec);
  }
#pragma omp barrier

  vec bigg1 = zeros<vec>(nr);
  vec bigg2 = zeros<vec>(nr);
  vec biggmat = zeros<vec>(nr);

#pragma omp parallel for
  for (int tp=1; tp < datapoints; tp++)
  {
    // set up a Gaussian pdf matrix as a column vector
    vec bigg1 = gaussian_pdf_vec_last(C1(tp), mu1, sig1);
    vec bigg2 = gaussian_pdf_vec_last(C2(tp), mu2, sig2);

    vec biggmat = bigg1 % bigg2;

    outpdf(tp-1) = k*k*dot(biggmat, approxpdfvec.col(tp-1));
  }
#pragma omp barrier 
  return outpdf;
}

// mat PDFcheck(const vec &thetavec, double h, double k, double yM)
// {
//   double h12 = sqrt(h);
//   int M = ceil(yM/k);
//   int veclen = 2*M+1;
//   int nr = veclen*veclen;
//   vec xvec = k*linspace<vec>(-M,M,veclen);

//   double supg = 0.5;
//   int gamma = ceil(5*h12*supg/k);
//   int dimg = 2*gamma;

//   mat Gnorm = zeros<mat>(veclen-dimg,veclen-dimg);

// #pragma omp parallel for
//   for (int newspace=0; newspace<nr; newspace++)
//   {
//     int j = floor(newspace/veclen); // subtract M to get math i, j \in [-M,M]
//     int i = newspace - veclen*floor(newspace/veclen);      
//     for (int ip=gamma; ip<(veclen-gamma); ip++)
//     {
//       for (int jp=gamma; jp<(veclen-gamma); jp++)
//       {
//         int test = jp*veclen + ip;
//         if ((test >= 0) && (test < nr))
//         {
//           {
//             double y1 = (ip-M)*k;
//             double y2 = (jp-M)*k;
//             double mu1 = y1 + f1(y1,y2,thetavec)*h;
//             double mu2 = y2 + f2(y1,y2,thetavec)*h;
//             double sigma1 = h12*g1(y1,y2,thetavec);
//             double sigma2 = h12*g2(y1,y2,thetavec);
//             double locbigg1 = gaussian_pdf(xvec(i), mu1, sigma1);
//             double locbigg2 = gaussian_pdf(xvec(j), mu2, sigma2);
//             Gnorm(ip-gamma,jp-gamma) += k*k*locbigg1*locbigg2;
//           }
//         }
//       }
//     }
//   }
// #pragma omp barrier
//   return(min(Gnorm));
// }

SEXP dtq2dCPP(SEXP s_thetavec, SEXP s_c1, SEXP s_c2, SEXP s_h, SEXP s_numsteps, SEXP s_k, SEXP s_yM)
{
    vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
    vec C1 = Rcpp::as<arma::vec>(s_c1);
    vec C2 = Rcpp::as<arma::vec>(s_c2);

    double h = Rcpp::as<double>(s_h);
    int numsteps = Rcpp::as<int>(s_numsteps);
    double k = Rcpp::as<double>(s_k);
    double yM = Rcpp::as<double>(s_yM);

    vec output = dtq(thetavec, C1, C2, h, numsteps, k, yM);
    return Rcpp::wrap( output );
}


// SEXP GCPP(SEXP s_thetavec, SEXP s_h, SEXP s_k, SEXP s_yM)
// {
//     vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
//     double h = Rcpp::as<double>(s_h);
//     double k = Rcpp::as<double>(s_k);
//     double yM = Rcpp::as<double>(s_yM);

//     vec vals = PDFcheck(thetavec, h, k, yM);
//     return Rcpp::wrap( vals );
// }

