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
    double p = diff_const * gaussian_pdf(x, mu, sigma);
    return p;
}

// gradient of the gaussian pdf, vectorized version of GSL function
static inline vec gradient1_gaussian_vec(const vec& x, const double mu, const double sigma)
{
    vec diff_const = (x - mu) / (2 * pow(fabs(sigma), 2));
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

// drift function
static inline double f(double x, const vec& thetavec)
{
  return(0.5*(thetavec(0) - x));
}

// diffusion function
static inline double g(double x, const vec& thetavec)
{
  return(thetavec(1));
}

mat gDTQ(const vec &thetavec, double h, double k, int M, int littlet, vec &init_data)
{
  // data is in init_data and is the vector of initial conditions for DTQ
  double h12 = sqrt(h);
  int numsteps = ceil(littlet/h);

  int veclen = 2*M + 1;
  int datapoints = init_data.n_elem;
  vec xvec = k*linspace<vec>(-M,M,veclen);

  mat A(veclen, datapoints);
  mat D1(veclen, datapoints);
  mat D2(veclen, datapoints);

  // compute the A and D matrix 
  for (int i = 0; i < datapoints; i++)
  {
    A.col(i) = gaussian_pdf_vec(xvec, init_data(i) + h*f(init_data(i), thetavec), h12*g(init_data(i), thetavec));
    D1.col(i) = gradient1_gaussian_vec(xvec, init_data(i) + h*f(init_data(i), thetavec), h12*g(init_data(i), thetavec));
    D2.col(i) = gradient2_gaussian_vec(xvec, init_data(i) + h*f(init_data(i), thetavec), h12*g(init_data(i), thetavec));
  }

  mat pdfmatrix(veclen, datapoints - 1);
  mat qmattheta1(veclen, datapoints - 1);
  mat qmattheta2(veclen, datapoints - 1);

  // loop over the timesteps
  for (int step = 1; step < numsteps - 1; step++)
  {
    mat Anew = zeros<mat>(veclen, datapoints);
    mat D1new = zeros<mat>(veclen, datapoints);
    mat D2new = zeros<mat>(veclen, datapoints);

    for (int i = 0; i < veclen; i++)
    {
      double y = xvec(i);

      double mu = y + h*f(y,thetavec);
      double sigma = h12*g(y,thetavec);

      double G = gaussian_pdf(xvec(i), mu, sigma);

      double Ggrad1 = gradient1_gaussian(xvec(i), mu, sigma);
      double Ggrad2 = gradient2_gaussian(xvec(i), mu, sigma);

      Anew.row(i) = k * G * A.row(i);
      D1new.row(i) = k * Ggrad1 * A.row(i) + k * G * D1.row(i);
      D2new.row(i) = k * Ggrad2 * A.row(i) + k * G * D2.row(i);
    }
    A = Anew;
    D1 = D1new;
    D2 = D2new;
  }
  return A;
}


SEXP gdtqCPP(SEXP s_thetavec, SEXP s_h, SEXP s_k, SEXP s_M, SEXP s_littlet, SEXP s_init_data)
{
  vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
  double h = Rcpp::as<double>(s_h);
  double k = Rcpp::as<double>(s_k);
  int M = Rcpp::as<int>(s_M);
  int littlet = Rcpp::as<int>(s_littlet);
  vec init_data = Rcpp::as<arma::vec>(s_init_data);

  mat mymat = gDTQ(thetavec, h, k, M, littlet, init_data);
  return Rcpp::wrap( mymat );
}


