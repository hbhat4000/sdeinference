#include "Rgdtq.h"

#include <RcppArmadillo.h>
#include <cmath>

using namespace std;
using namespace arma;

// drift function
static inline double f(double x, const vec& thetavec)
{
  return(thetavec(0)*(x)*(thetavec(1) - pow(x, 2)));
}

// diffusion function
static inline double g(double x, const vec& thetavec)
{
  return(thetavec(2));
}

// Function overloading
static inline double gaussian_pdf(int choice, const double y, const double h, const vec& thetavec, double driftval, double diffval)
{
  double h12 = sqrt(h);

  double mu = y + h * driftval;
  double sigma = fabs(h12 * diffval);
  double part1 = (y - mu);
  double u = part1 / sigma;
  double G = (1 / (sqrt (2 * M_PI) * sigma)) * exp (-(u * u) / 2);

  double dGdf = G * part1 / (pow(diffval, 2));
  double dGdg = G * (-pow(diffval, -1) + pow(part1, 2) / (pow(diffval, 3) * h));

  if(choice == 0)
  {
    return G;
  }
  else if(choice == 1)
  {
    double dfdtheta = driftval / thetavec(0);
    double dgdtheta = 0;
    double grad = dGdf * dfdtheta + dGdg * dgdtheta;
    return grad;
  }
  else if(choice == 2)
  {
    double dfdtheta = thetavec(0) * y;
    double dgdtheta = 0;
    double grad = dGdf * dfdtheta + dGdg * dgdtheta;
    return grad;
  }
  else if(choice == 3)
  {
    double dfdtheta = 0;
    double dgdtheta = 1;
    double grad = dGdf * dfdtheta + dGdg * dgdtheta;
    return grad;
  }
  else
  {
    cout << "Incorrect option" << endl;
    return 0;
  }
}

static inline vec gaussian_pdf(int choice, const vec& x, const double y, const double h, const vec& thetavec, double driftval, double diffval)
{
  double h12 = sqrt(h);

  double mu = y + h * driftval;
  double sigma = fabs(h12 * diffval);
  vec part1 = (x - mu);
  vec u = part1 / sigma;
  vec G = (1 / (sqrt (2 * M_PI) * sigma)) * exp (-(u % u) / 2);

  vec dGdf = (G % part1) / (pow(diffval, 2));
  vec dGdg = G % (-1/diffval + (part1 % part1) / (pow(diffval, 3) * h));

  if(choice == 0)
  {
    return G;
  }
  else if(choice == 1)
  {
    double dfdtheta = driftval / thetavec(0);
    double dgdtheta = 0;
    vec grad = dGdf * dfdtheta + dGdg * dgdtheta;
   return grad;
  }
  else if(choice == 2)
  {
    double dfdtheta = thetavec(0) * y;
    double dgdtheta = 0;
    vec grad = dGdf * dfdtheta + dGdg * dgdtheta;
    return grad;
  }
  else if(choice == 3)
  {
    double dfdtheta = 0;
    double dgdtheta = 1;
    vec grad = dGdf * dfdtheta + dGdg * dgdtheta;
    return grad;
  }
  else
  {
    cout << "Incorrect option" << endl;
    return 0;
  }
}

// Ornstein Uhlenbeck SDE
// dX_t = theta1 X_t (theta2 - (X_t)^2) dt +  theta3 dW_t
// dx_t = f dt + g dW_t

mat gDTQ(const vec &thetavec, double h, double k, int M, int littlet, vec &init_data)
{
  int numsteps = ceil(littlet/h);

  int veclen = 2*M + 1;
  int datapoints = init_data.n_elem;
  vec xvec = k*linspace<vec>(-M,M,veclen);

  mat A(veclen, datapoints);
  mat D1(veclen, datapoints);
  mat D2(veclen, datapoints);
  mat D3(veclen, datapoints);

  int choice = 0;

  // compute the A and D matrix 
  for (int i = 0; i < datapoints; i++)
  {
    double y = init_data(i);
    double driftval = f(y, thetavec);
    double diffval = g(y, thetavec);

    A.col(i) = gaussian_pdf(choice = 0, xvec, y, h, thetavec, driftval, diffval);
    D1.col(i) = gaussian_pdf(choice = 1, xvec, y, h, thetavec, driftval, diffval);
    D2.col(i) = gaussian_pdf(choice = 2, xvec, y, h, thetavec, driftval, diffval);
    D3.col(i) = gaussian_pdf(choice = 3, xvec, y, h, thetavec, driftval, diffval);
    // A.col(i) = gaussian_pdf_vec(xvec, init_data(i) + h*f(init_data(i), thetavec), h12*g(init_data(i), thetavec));
    // D1.col(i) = gradient1_gaussian_vec(xvec, init_data(i) + h*f(init_data(i), thetavec), h12*g(init_data(i), thetavec));
    // D2.col(i) = gradient2_gaussian_vec(xvec, init_data(i) + h*f(init_data(i), thetavec), h12*g(init_data(i), thetavec));
  }

  mat pdfmatrix(veclen, datapoints - 1);
  mat qmattheta1(veclen, datapoints - 1);
  mat qmattheta2(veclen, datapoints - 1);
  mat qmattheta3(veclen, datapoints - 1);

  // loop over the timesteps
  for (int step = 1; step < numsteps - 1; step++)
  {
    mat Anew = zeros<mat>(veclen, datapoints);
    mat D1new = zeros<mat>(veclen, datapoints);
    mat D2new = zeros<mat>(veclen, datapoints);
    mat D3new = zeros<mat>(veclen, datapoints);

    for (int i = 0; i < veclen; i++)
    {
      double y = xvec(i);
      double driftval = f(y, thetavec);
      double diffval = g(y, thetavec);

      double G = gaussian_pdf(choice = 0, y, h, thetavec, driftval, diffval);
      double Ggrad1 = gaussian_pdf(choice = 1, y, h, thetavec, driftval, diffval);
      double Ggrad2 = gaussian_pdf(choice = 2, y, h, thetavec, driftval, diffval);
      double Ggrad3 = gaussian_pdf(choice = 3, y, h, thetavec, driftval, diffval);
      // double G = gaussian_pdf(xvec(i), mu, sigma);

      // double Ggrad1 = gradient1_gaussian(xvec(i), mu, sigma);
      // double Ggrad2 = gradient2_gaussian(xvec(i), mu, sigma);

      Anew.row(i) = k * G * A.row(i);
      D1new.row(i) = k * Ggrad1 * A.row(i) + k * G * D1.row(i);
      D2new.row(i) = k * Ggrad2 * A.row(i) + k * G * D2.row(i);
      D3new.row(i) = k * Ggrad3 * A.row(i) + k * G * D3.row(i);
    }
    A = Anew;
    D1 = D1new;
    D2 = D2new;
    D3 = D3new;
  }

  mat temp = join_cols(D1, D2);
  mat temp2 = join_cols(temp, D3);
  mat returningmat = join_cols(A, temp2);
  return returningmat;
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


// // gaussian pdf, non-vectorized version 
// static inline double gaussian_pdf(const double x, const double mu, const double sigma)
// {
//   double u = (x - mu) / fabs(sigma);
//   double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u * u) / 2);
//   return p;
// }

// // gaussian pdf, vectorized version 
// static inline vec gaussian_pdf_vec(const vec& x, const double mu, const double sigma)
// {
//   vec u = (x - mu) / fabs(sigma);
//   vec p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u % u) / 2);
//   return p;
// }

// // gradient of the gaussian pdf, non-vectorized version 
// static inline double gradient1_gaussian(const double x, const double mu, const double sigma)
// {
//   double diff_const = (x - mu) / (2 * pow(fabs(sigma), 2));
//   double p = diff_const * gaussian_pdf(x, mu, sigma);
//   return p;
// }

// // gradient of the gaussian pdf, vectorized version 
// static inline vec gradient1_gaussian_vec(const vec& x, const double mu, const double sigma)
// {
//   vec diff_const = (x - mu) / (2 * pow(fabs(sigma), 2));
//   vec p = diff_const % gaussian_pdf_vec(x, mu, sigma);
//   return p;
// }

// // gradient of the gaussian pdf, non-vectorized version
// static inline double gradient2_gaussian(const double x, const double mu, const double sigma)
// {
//   double diff_const = (pow((x - mu), 2) - pow(fabs(sigma), 2)) / pow(fabs(sigma),3);
//   double p = diff_const * gaussian_pdf(x, mu, sigma);
//   return p;
// }

// // gradient of the gaussian pdf, vectorized version 
// static inline vec gradient2_gaussian_vec(const vec& x, const double mu, const double sigma)
// {
//   vec diff_const = ((x - mu) % (x - mu) - pow(fabs(sigma), 2)) / pow(fabs(sigma),3);
//   vec p = diff_const % gaussian_pdf_vec(x, mu, sigma);
//   return p;
// }

// // gradient of the gaussian pdf, non vectorized version
// static inline double gradient3_gaussian(const double x, const double mu, const double sigma)
// {
//   double diff_const = (pow((x - mu), 2) - pow(fabs(sigma), 2)) / pow(fabs(sigma),3);
//   double p = diff_const * gaussian_pdf(x, mu, sigma);
//   return p;
// }

// // gradient of the gaussian pdf, vectorized version
// static inline vec gradient3_gaussian_vec(const vec& x, const double mu, const double sigma)
// {
//   vec diff_const = ((x - mu) % (x - mu) - pow(fabs(sigma), 2)) / pow(fabs(sigma),3);
//   vec p = diff_const % gaussian_pdf_vec(x, mu, sigma);
//   return p;
// }
