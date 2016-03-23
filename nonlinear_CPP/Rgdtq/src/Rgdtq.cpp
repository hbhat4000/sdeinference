#include "Rgdtq.h"

#include <RcppArmadillo.h>
#include <cmath>

using namespace std;
using namespace arma;

// gaussian pdf, non-vectorized version of GSL function
static inline double gaussian_pdf(const double x, const double mu, const double sigma)
{
    double u = (x - mu) / fabs(sigma);
    double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u * u) / 2);
    return p;
}

// gaussian pdf, vectorized version of GSL function
static inline vec gaussian_pdf_vec(const vec& x, const double mu, const double sigma)
{
    vec u = (x - mu) / fabs(sigma);
    vec p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u % u) / 2);
    return p;
}

// gradient of the gaussian pdf, non-vectorized version of GSL function
static inline double gradient1_gaussian(const double x, const double mu, const double sigma)
{
    double diff_const = (x - mu) / (2 * pow(fabs(sigma), 2));
    // double p = (diff_const / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u * u) / 2);
    double p = diff_const * gaussian_pdf(x, mu, sigma);
    return p;
}

// gradient of the gaussian pdf, vectorized version of GSL function
static inline vec gradient1_gaussian_vec(const vec& x, const double mu, const double sigma)
{
    vec diff_const = (x - mu) / (2 * pow(fabs(sigma), 2));
    // vec diff_const = u / (2 * fabs(sigma));
    vec p = diff_const % gaussian_pdf_vec(x, mu, sigma);
    return p;
}

// gradient of the gaussian pdf, non-vectorized version of GSL function
static inline double gradient2_gaussian(const double x, const double mu, const double sigma)
{
    double diff_const = (pow((x - mu), 2) - pow(fabs(sigma), 2)) / pow(fabs(sigma),3);
    double p = diff_const * gaussian_pdf(x, mu, sigma);
    return p;
}

// gradient of the gaussian pdf, vectorized version of GSL function
static inline vec gradient2_gaussian_vec(const vec& x, const double mu, const double sigma)
{
    vec diff_const = ((x - mu) % (x - mu) - pow(fabs(sigma), 2)) / pow(fabs(sigma),3);
    vec p = diff_const % gaussian_pdf_vec(x, mu, sigma);
    return p;
}

// Ornstein Uhlenbeck SDE
// dx_t = 0.5*(theta_1 - x_t)dt + theta_2 dW_t
// dx_t = f dt + g dW_t

// f function
static inline double f(double x, const vec& thetavec)
{
  return(0.5*(thetavec(0) - x));
}

// g function
static inline double g(double x, const vec& thetavec)
{
  return(thetavec(1));
}

mat gDTQ(const vec &thetavec, const vec &C0, double h, int numsteps, double k, double yM)
{
  // each element in C0 represents an initial condition of the DTQ method
  double h12 = sqrt(h);
  int M = ceil(yM/k);
  int veclen = 2*M+1;
  int datapoints = C0.n_elem;
  vec xvec = k*linspace<vec>(-M,M,veclen);
  // mat Phat(veclen, datapoints);
  // mat Pgrad1(veclen, datapoints);
  // mat Pgrad2(veclen, datapoints);
  mat P(3*veclen, datapoints);

  for (int i = 0; i < datapoints; i++)
  {
    // Phat.col(i) = gaussian_pdf_vec(xvec, C0(i) + h*f(C0(i), thetavec), h12*g(C0(i),thetavec));
    // Pgrad1.col(i) = gradient1_gaussian_vec(xvec, C0(i) + h*f(C0(i), thetavec), h12*g(C0(i),thetavec));
    // Pgrad2.col(i) = gradient2_gaussian_vec(xvec, C0(i) + h*f(C0(i), thetavec), h12*g(C0(i),thetavec));
  }

  // double supg = max(g(0,thetavec));
  // int gamma = ceil(5*h12*supg/k);

  // loop over the timesteps
  for (int step = 1; step < numsteps - 1; step++)
  {
    // mat Phatnew = zeros<mat>(veclen, datapoints);
    // mat Pgrad1new = zeros<mat>(veclen, datapoints);
    // mat Pgrad2new = zeros<mat>(veclen, datapoints);
    mat Pnew = zeros<mat>(3*veclen, datapoints);
    // Index 0:veclen-1 Phat
    // Index veclen:2veclen - 1 Pgrad1
    // Index 2veclen:3veclen -1 Pgrad2

    for (int i = 0; i < veclen; i++)
    {
      // only do the computation gamma away from boundary
      double y = xvec(i);
      double mu = y + h*f(y,thetavec);
      double sigma = h12*g(y,thetavec);
      double G = gaussian_pdf(xvec(i), mu, sigma);
      double Ggrad1 = gradient1_gaussian(xvec(i), mu, sigma);
      double Ggrad2 = gradient2_gaussian(xvec(i), mu, sigma);
      // Phatnew.row(i) = k*G*Phat.row(i);
      // Pgrad1new.row(i) = k*Ggrad1*Phat.row(i) + k*G*Pgrad1.row(i);
      // Pgrad2new.row(i) = k*Ggrad2*Phat.row(i) + k*G*Pgrad2.row(i);

      Pnew.row(i) = k*G*P.row(i);
      Pnew.row(veclen + i) = k*Ggrad1*P.row(i) + k*G*P.row(veclen + i);
      Pnew.row(2*veclen + i) = k*Ggrad2*P.row(i) + k*G*P.row(2*veclen + i);
    }
    // Phat = Phatnew;
    // Pgrad1 = Pgrad1new;
    // Pgrad2 = Pgrad2new;
    P = Pnew
  }
  return P;
}


SEXP gdtqCPP(SEXP s_thetavec, SEXP s_c0, SEXP s_h, SEXP s_numsteps, SEXP s_k, SEXP s_yM)
{
    vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
    vec C0 = Rcpp::as<arma::vec>(s_c0);
    double h = Rcpp::as<double>(s_h);
    int numsteps = Rcpp::as<int>(s_numsteps);
    double k = Rcpp::as<double>(s_k);
    double yM = Rcpp::as<double>(s_yM);

    mat mymat = gDTQ(thetavec, C0, h, numsteps, k, yM);
    return Rcpp::wrap( mymat );
}

// SEXP GCPP(SEXP s_thetavec, SEXP s_h, SEXP s_k, SEXP s_yM)
// {
//     vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
//     double h = Rcpp::as<double>(s_h);
//     double k = Rcpp::as<double>(s_k);
//     double yM = Rcpp::as<double>(s_yM);

//     double Gval = PDFcheck(thetavec, h, k, yM);
//     return Rcpp::wrap( Gval );
// }


