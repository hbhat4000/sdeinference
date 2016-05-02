#include "Rgdtq.h"

#include <RcppArmadillo.h>
#include <cmath>

using namespace std;
using namespace arma;

/* Linear case */
// Ornstein - Uhlenback (linear case)
static inline mat f(const mat& Y, const vec& thetavec)
{
  return(thetavec(0) * (thetavec(1) - Y));
}

static inline mat g(const mat& Y, const vec& thetavec)
{
  mat temp = Y;
  return(temp.fill(thetavec(2)));
}

static inline cube Dtheta(const vec &xvec, const vec &yvec, const double h, const vec &thetavec)
{
  int ny = yvec.n_elem;
  int nx = xvec.n_elem;
  mat X(nx, ny);
  mat temp(ny, nx);

  X.each_col() = xvec;
  temp.each_col() = yvec;
  mat Y = trans(temp);
  mat thisdrift = f(Y, thetavec);
  mat thisdiff = abs(g(Y, thetavec));
  mat sqthisdiff = thisdiff % thisdiff;
  mat part1 = X - Y - thisdrift*h;

  mat G = exp(- part1 % part1 / (2 * h * sqthisdiff)) / (sqrt(2 * M_PI * h * sqthisdiff));

  mat dGdf = (G % part1) / sqthisdiff;
  mat dGdg = G % (-1/thisdiff + (part1 % part1) / (sqthisdiff % thisdiff * h));

  cube Ggrad(nx, ny, thetavec.n_elem + 1);

  Ggrad.slice(0) = G;
  Ggrad.slice(1) = dGdf % (thisdrift / thetavec(0));
  Ggrad.slice(2) = dGdf * thetavec(0);
  Ggrad.slice(3) = dGdg;

  return Ggrad;
}

cube gDTQ(const vec &thetavec, const double h, const double k, const int M, const int littlet, const mat &init_data)
{
  int numsteps = ceil(littlet/h);

  int veclen = 2*M + 1;
  int datapoints = init_data.n_cols;
  vec xvec = k*(linspace<vec>(-M,M,veclen));
  int numtheta = thetavec.n_elem;
  
// D.slice(0) = A, D.slice(1) = Dtheta1, D.slice(2) = Dtheta2, D.slice(3) = Dtheta3
// qmattheta.slice(0) = pdfmatrix, rest qmattheta.slice(i)

  cube D = Dtheta(xvec, xvec, h, thetavec);
  cout << "Dimensions of cube D: " << D.n_rows << ", " << D.n_cols << ", " << D.n_slices << endl;
  cout << "Sum of cube D: " << accu(D) << endl;
  // cout << "3*3 matrix for objective function: " << D.subcube(0, 0, 0, 2, 2, 3) << endl;

  cube qmattheta = zeros<cube>(veclen, datapoints - 1, numtheta + 1);
  
  for(int curcol = 0; curcol < datapoints - 1; curcol++)
  {
    cube initderivs = Dtheta(xvec, init_data.col(curcol), h, thetavec);
    if(curcol == 10)
    {
      cout << "Dimensions of cube initderivs:" << initderivs.n_rows << ", " << initderivs.n_cols << ", " << initderivs.n_slices << endl;
      cout << "Sum of cube initderivs: " << accu(initderivs) << endl;
    }

    for(int i = 0; i < numtheta + 1; i++)
    {
      (qmattheta.slice(i)).col(curcol) = mean(initderivs.slice(i), 1);
    }
  }
  int startstep = 1;
  for(int curstep = startstep; curstep < numsteps - 1; curstep++)
  {
    qmattheta.slice(0) = k * D.slice(0) * qmattheta.slice(0);

    for(int i = 1; i < numtheta + 1; i++)
    {
      qmattheta.slice(i) = k * D.slice(0) * qmattheta.slice(i) + k * D.slice(i) * qmattheta.slice(0);
    }
  }
  cout << "Dimensions of cube qmattheta:" << qmattheta.n_rows << ", " << qmattheta.n_cols << ", " << qmattheta.n_slices << endl;
  cout << "Sum of cube qmattheta: " << accu(qmattheta) << endl;

// gradient.slice(0) = likelihood
// gdmat.slice(0) = gammamat
  cube gradient = zeros<cube>(init_data.n_rows, datapoints - 1, numtheta + 1);
  for(int curcol = 1; curcol < datapoints; curcol++)
  {
    cube gdmat = Dtheta(init_data.col(curcol), xvec, h, thetavec);
    // if(curcol == 1)
    {
      cout << "Dimensions of cube gdmat:" << gdmat.n_rows << ", " << gdmat.n_cols << ", " << gdmat.n_slices << endl;
      cout << "Sum of cube gdmat: " << accu(gdmat) << endl;
    // cout << "gdmat for curcol " << curcol << ": " << gdmat << endl;
    }

    (gradient.slice(0)).col(curcol - 1) = k * gdmat.slice(0) * (qmattheta.slice(0)).col(curcol - 1);

    for(int i = 1; i < numtheta + 1; i++)
    {
  cout << "Dimensions of cube gradient:" << gradient.n_rows << ", " << gradient.n_cols << ", " << gradient.n_slices << endl;
  cout << "Sum of cube gradient: " << accu(gradient) << endl;
      (gradient.slice(i)).col(curcol - 1) = k * gdmat.slice(i) * (qmattheta.slice(0)).col(curcol - 1) + k * (gdmat.slice(0)) * (qmattheta.slice(i)).col(curcol - 1);
    }
  }
  cout << "Dimensions of cube gradient:" << gradient.n_rows << ", " << gradient.n_cols << ", " << gradient.n_slices << endl;
  cout << "Sum of cube gradient: " << accu(gradient) << endl;

  // // to check that the right values are getting passed to R
  cout << "objective: " << -accu(log(gradient.slice(0))) << endl;
  cout << "gradient 1: " << -accu(gradient.slice(1)/gradient.slice(0)) << endl;
  cout << "gradient 2: " << -accu(gradient.slice(2)/gradient.slice(0)) << endl;
  cout << "gradient 3: " << -accu(gradient.slice(3)/gradient.slice(0)) << endl;


return gradient;
}

static inline mat Dthetadtq(const vec &xvec, const vec &yvec, const double h, const vec &thetavec)
{
  int ny = yvec.n_elem;
  int nx = xvec.n_elem;
  mat X(nx, ny);
  mat temp(ny, nx);

  X.each_col() = xvec;
  temp.each_col() = yvec;
  mat Y = trans(temp);
  mat thisdrift = f(Y, thetavec);
  mat thisdiff = abs(g(Y, thetavec));
  mat sqthisdiff = 2 * h * square(thisdiff);
  mat part1 = X - Y - thisdrift*h;

  mat G = exp(-square(part1) / sqthisdiff) / (sqrt(M_PI * sqthisdiff));
  return G;
}

mat DTQ(const vec &thetavec, const double h, const double k, const int M, const int littlet, const mat &init_data)
{
  int numsteps = ceil(littlet/h);
  int veclen = 2*M + 1;
  int datapoints = init_data.n_cols;
  vec xvec = k*linspace<vec>(-M, M, veclen);

  mat D = Dthetadtq(xvec, xvec, h, thetavec);
  mat pdfmatrix = zeros<mat>(veclen, datapoints - 1);
  
  for(int curcol = 0; curcol < datapoints - 1; curcol++)
  {
    mat initderivs = Dthetadtq(xvec, init_data.col(curcol), h, thetavec);
    pdfmatrix.col(curcol) = mean(initderivs, 1);       // mean(,1) for rowMeans
  }

  int startstep = 1;
  for(int curstep = startstep; curstep < numsteps - 1; curstep++)
  {
    pdfmatrix = k * D * pdfmatrix;
  }

  mat likelihood = zeros<mat>(init_data.n_rows, datapoints - 1);
  for(int curcol = 1; curcol < datapoints; curcol++)
  {
    mat gdmat = Dthetadtq(init_data.col(curcol), xvec, h, thetavec);
    likelihood.col(curcol - 1) = k * gdmat * pdfmatrix.col(curcol - 1);
  }

return likelihood;
}

SEXP gdtqCPP_linear(SEXP s_thetavec, SEXP s_h, SEXP s_k, SEXP s_M, SEXP s_littlet, SEXP s_init_data)
{
  vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
  double h = Rcpp::as<double>(s_h);
  double k = Rcpp::as<double>(s_k);
  int M = Rcpp::as<int>(s_M);
  int littlet = Rcpp::as<int>(s_littlet);
  mat init_data = Rcpp::as<arma::mat>(s_init_data);

  cube grad = gDTQ(thetavec, h, k, M, littlet, init_data);
  return Rcpp::wrap( grad );
}

SEXP dtqCPP_linear(SEXP s_thetavec, SEXP s_h, SEXP s_k, SEXP s_M, SEXP s_littlet, SEXP s_init_data)
{
  vec thetavec = Rcpp::as<arma::vec>(s_thetavec);
  double h = Rcpp::as<double>(s_h);
  double k = Rcpp::as<double>(s_k);
  int M = Rcpp::as<int>(s_M);
  int littlet = Rcpp::as<int>(s_littlet);
  mat init_data = Rcpp::as<arma::mat>(s_init_data);

  mat lik = DTQ(thetavec, h, k, M, littlet, init_data);
  return Rcpp::wrap( lik );
}
