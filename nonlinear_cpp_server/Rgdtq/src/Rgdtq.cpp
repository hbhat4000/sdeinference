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

static inline vec f(vec x, const vec& thetavec)
{
  return(thetavec(0)*(x)*(thetavec(1) - pow(x, 2)));
}

// diffusion function
static inline double g(double x, const vec& thetavec)
{
  return(thetavec(2));
}

static inline mat integrandmat(const vec &xvec, const vec &yvec, const double h, const vec &thetavec)
{
  int ny = yvec.n_elem;
  int nx = xvec.n_elem;
  mat X(nx, ny);
  mat temp(ny, nx);
  mat Y;
  mat thisdrift;
  Y.copy_size(X);
  thisdrift.copy_size(X);

  X.each_col() = xvec;
  temp.each_col() = yvec;
  Y = temp.t();

  double thisdiff = abs(thetavec(2));
  thisdrift = thetavec(0) * Y % (thetavec(1) - Y % Y);

  mat out(nx, ny);
  out = exp(- (X - Y - thisdrift*h) % (X - Y - thisdrift*h) / (2 * h * (thisdiff * thisdiff))) / (thisdiff * sqrt(2 * M_PI * h));
  return out;
}

static inline cube Dtheta(const vec &xvec, const vec &yvec, const double h, const vec &thetavec)
{
  int ny = yvec.n_elem;
  int nx = xvec.n_elem;
  mat G = integrandmat(xvec, yvec, h, thetavec);
  mat X(nx, ny);
  mat temp(ny, nx);

  X.each_col() = xvec;
  temp.each_col() = yvec;
  mat Y = temp.t();
  mat thisdrift = thetavec(0) * Y % (thetavec(1) - Y % Y);
  double thisdiff = abs(thetavec(2));
  mat part1 = X - Y - thisdrift*h;
  mat dGdf = (G % part1) / pow(thisdiff, 2);
  mat dGdg = G % (-1/thisdiff + (part1 % part1) / (pow(thisdiff, 3) * h));
  cube grad(nx, ny, thetavec.n_elem);
  grad.slice(0) = dGdf % (thisdrift / thetavec(0));
  grad.slice(1) = dGdf % (thetavec(0) * Y);
  grad.slice(2) = dGdg;

  return grad;
}

cube gDTQ(const vec &thetavec, double h, double k, int M, int littlet, mat &init_data)
{
  int numsteps = ceil(littlet/h);

  int veclen = 2*M + 1;
  int datapoints = init_data.n_cols;
  vec xvec = k*linspace<vec>(-M,M,veclen);
  int numtheta = thetavec.n_elem;

  mat A = integrandmat(xvec, xvec, h, thetavec);
  cube D = Dtheta(xvec, xvec, h, thetavec);
  mat pdfmatrix(veclen, datapoints - 1);

  cube qmattheta = zeros<cube>(veclen, datapoints - 1, numtheta);
  cube initderivs;
  initderivs.copy_size(qmattheta);

  for(int curcol = 0; curcol < datapoints - 1; curcol++)
  {
    mat B = integrandmat(xvec, init_data.col(curcol), h, thetavec);

    pdfmatrix.col(curcol) = sum(B, 1);

    initderivs = Dtheta(xvec, init_data.col(curcol), h, thetavec);

    for(int i = 0; i < numtheta; i++)
    {
      (qmattheta.slice(i)).col(curcol) = sum(initderivs.slice(i), 1);
    }
  }

  int startstep = 1;
  for(int curstep = startstep; curstep < numsteps - 1; curstep++)
  {
    pdfmatrix = k * A * pdfmatrix;
    for(uword i = 0; i < qmattheta.n_slices; i++)
    {
      qmattheta.slice(i) = k * A * qmattheta.slice(i) + k * D.slice(i) * pdfmatrix;
    }
  }

  mat likelihood = zeros<mat>(init_data.n_rows, datapoints - 1);
  cube gradient = zeros<cube>(init_data.n_rows, datapoints - 1, numtheta);

  for(int curcol = 1; curcol < datapoints; curcol++)
  {
    mat gammamat = integrandmat(init_data.col(curcol), xvec, h, thetavec);
    likelihood.col(curcol - 1) = k * gammamat * pdfmatrix.col(curcol - 1);
    cube gdmat = Dtheta(init_data.col(curcol), xvec, h, thetavec);

    for(uword i = 0; i < gradient.n_slices; i++)
    {
      (gradient.slice(i)).col(curcol - 1) = k * gdmat.slice(i) * pdfmatrix.col(curcol - 1) + k * gammamat * (qmattheta.slice(i)).col(curcol - 1);
    }
  }

  cube out(init_data.n_rows, datapoints - 1, numtheta + 1);
  out.slice(0) = likelihood;
  for(int i = 1; i < numtheta + 1; i++)
  {
    out.slice(i) = gradient.slice(i-1);
  }

return out;
}
  // int choice;
  // cube pdf = zeros<cube>(veclen, datapoints - 1, numtheta + 1);
  // cube pdfnew = zeros<cube>(veclen, datapoints - 1, numtheta + 1);
  // compute the A and D matrix 
  // for(int j = 0; j < numtheta + 1; j++)
  // {
  //   for (int i = 0; i < datapoints; i++)
  //   {
  //     double y = init_data(i);
  //     double driftval = f(y, thetavec);
  //     double diffval = g(y, thetavec);

  //     (pdf.slice(j)).col(i) = gaussian_pdf(choice = j, xvec, y, h, thetavec, driftval, diffval);
  //   }
  // }
  // cout << "output 3" << endl;

//   // loop over the timesteps
//   for (int step = 1; step < numsteps - 1; step++)
//   {
//     for (int i = 0; i < datapoints; i++)
//     {
//       double y = xvec(i);
//       double driftval = f(y, thetavec);
//       double diffval = g(y, thetavec);

//       vec Ggrad(numtheta + 1);

//       for(int j = 0; j < numtheta + 1; j++)
//       {
//         Ggrad(j) = gaussian_pdf(choice = j, y, h, thetavec, driftval, diffval);
//         if(j == 0)
//           (pdfnew.slice(0)).col(i) = k * Ggrad(0) * (pdf.slice(0)).col(i);
//         else
//           (pdfnew.slice(j)).col(i) = k * Ggrad(j) * (pdf.slice(0)).col(i) + k * Ggrad(0) * (pdf.slice(j)).col(i);
//       }
//     }

//   pdf = pdfnew;
//   }
//   cout << "output 4" << endl;

//   return pdf;
// }


SEXP gdtqCPP(SEXP s_thetavec, SEXP s_h, SEXP s_k, SEXP s_M, SEXP s_littlet, SEXP s_init_data)
{
  vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
  double h = Rcpp::as<double>(s_h);
  double k = Rcpp::as<double>(s_k);
  int M = Rcpp::as<int>(s_M);
  int littlet = Rcpp::as<int>(s_littlet);
  mat init_data = Rcpp::as<arma::mat>(s_init_data);

  cube mymat = gDTQ(thetavec, h, k, M, littlet, init_data);
  return Rcpp::wrap( mymat );
}


// Function overloading
// static inline double gaussian_pdf(int choice, const double y, const double h, const vec& thetavec, double driftval, double diffval)
// {
//   double h12 = sqrt(h);

//   double mu = y + h * driftval;
//   double sigma = fabs(h12 * diffval);
//   double part1 = (y - mu);
//   double u = part1 / sigma;
//   double G = (1 / (sqrt (2 * M_PI) * sigma)) * exp (-(u * u) / 2);

//   double dGdf = G * part1 / (pow(diffval, 2));
//   double dGdg = G * (-pow(diffval, -1) + pow(part1, 2) / (pow(diffval, 3) * h));

//   if(choice == 0)
//   {
//     return G;
//   }
//   else if(choice == 1)
//   {
//     double dfdtheta = driftval / thetavec(0);
//     double dgdtheta = 0;
//     double grad = dGdf * dfdtheta + dGdg * dgdtheta;
//     return grad;
//   }
//   else if(choice == 2)
//   {
//     double dfdtheta = thetavec(0) * y;
//     double dgdtheta = 0;
//     double grad = dGdf * dfdtheta + dGdg * dgdtheta;
//     return grad;
//   }
//   else if(choice == 3)
//   {
//     double dfdtheta = 0;
//     double dgdtheta = 1;
//     double grad = dGdf * dfdtheta + dGdg * dgdtheta;
//     return grad;
//   }
//   else
//   {
//     cout << "Incorrect option" << endl;
//     return 0;
//   }
// }

// static inline vec gaussian_pdf(int choice, const vec& x, const double y, const double h, const vec& thetavec, double driftval, double diffval)
// {
//   double h12 = sqrt(h);

//   double mu = y + h * driftval;
//   double sigma = fabs(h12 * diffval);
//   vec part1 = (x - mu);
//   vec u = part1 / sigma;
//   vec G = (1 / (sqrt (2 * M_PI) * sigma)) * exp (-(u % u) / 2);

//   vec dGdf = (G % part1) / (pow(diffval, 2));
//   vec dGdg = G % (-1/diffval + (part1 % part1) / (pow(diffval, 3) * h));

//   if(choice == 0)
//   {
//     return G;
//   }
//   else if(choice == 1)
//   {
//     double dfdtheta = driftval / thetavec(0);
//     double dgdtheta = 0;
//     vec grad = dGdf * dfdtheta + dGdg * dgdtheta;
//    return grad;
//   }
//   else if(choice == 2)
//   {
//     double dfdtheta = thetavec(0) * y;
//     double dgdtheta = 0;
//     vec grad = dGdf * dfdtheta + dGdg * dgdtheta;
//     return grad;
//   }
//   else if(choice == 3)
//   {
//     double dfdtheta = 0;
//     double dgdtheta = 1;
//     vec grad = dGdf * dfdtheta + dGdg * dgdtheta;
//     return grad;
//   }
//   else
//   {
//     cout << "Incorrect option" << endl;
//     return 0;
//   }
// }



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
