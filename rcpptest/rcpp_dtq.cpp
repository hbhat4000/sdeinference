#include <Rcpp.h>
#include <cmath>
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

// [[Rcpp::export]]
cube Dtheta(NumericVector xvec, NumericVector yvec, double h, NumericVector thetavec)
{
  int ny = yvec.n_elem;
  int nx = xvec.n_elem;
  mat X(nx, ny);
  mat temp(ny, nx);

  X.each_col() = xvec;
  temp.each_col() = yvec;
  mat Y = trans(temp);
  mat thisdrift = thetavec(0) * (thetavec(1) - Y);
  mat temp = Y;
  mat thisdiff = abs(temp.fill(thetavec(2)));
  mat sqthisdiff = square(thisdiff);
  mat part1 = X - Y - thisdrift*h;

  mat G = exp(-square(part1) / (2 * h * sqthisdiff)) / (sqrt(2 * M_PI * h * sqthisdiff));

  mat dGdf = (G % part1) / sqthisdiff;
  mat dGdg = G % (-1/thisdiff + square(part1) / (sqthisdiff % thisdiff * h));

  cube Ggrad(nx, ny, thetavec.n_elem + 1);

  Ggrad.slice(0) = G;
  Ggrad.slice(1) = dGdf % (thisdrift / thetavec(0));
  Ggrad.slice(2) = dGdf * thetavec(0);
  Ggrad.slice(3) = dGdg;

  return Ggrad;
}

// [[Rcpp::export]]
NumericCube gDTQ(NumericVector thetavec, double h, double k, int M, int littlet, NumericMatrix init_data)
{
  int numsteps = ceil(littlet/h);

  int veclen = 2*M + 1;
  int datapoints = init_data.n_cols;
  vec xvec = k*linspace<vec>(-M,M,veclen);
  int numtheta = thetavec.n_elem;

// D.slice(0) = A, D.slice(1) = Dtheta1, D.slice(2) = Dtheta2, D.slice(3) = Dtheta3
// qmattheta.slice(0) = pdfmatrix, rest qmattheta.slice(i)

  cube D = Dtheta(xvec, xvec, h, thetavec);
  cube qmattheta = zeros<cube>(veclen, datapoints - 1, numtheta + 1);
  
  for(int curcol = 0; curcol < datapoints - 1; curcol++)
  {
    cube initderivs = Dtheta(xvec, init_data.col(curcol), h, thetavec);
    for(int i = 0; i < numtheta + 1; i++)
    {
      (qmattheta.slice(i)).col(curcol) = sum(initderivs.slice(i), 1);
    }
  }

  int startstep = 1;
  for(int curstep = startstep; curstep < numsteps - 1; curstep++)
  {
    for(int i = 0; i < numtheta + 1; i++)
    {
      if(i == 0)
      {
        qmattheta.slice(0) = k * D.slice(0) * qmattheta.slice(0);
      }
      else
      {
        qmattheta.slice(i) = k * D.slice(0) * qmattheta.slice(i) + k * D.slice(i) * qmattheta.slice(0);
      }
    }
  }

// gradient.slice(0) = likelihood
// gdmat.slice(0) = gammamat
  cube gradient = zeros<cube>(init_data.n_rows, datapoints - 1, numtheta + 1);
  for(int curcol = 1; curcol < datapoints; curcol++)
  {
    cube gdmat = Dtheta(init_data.col(curcol), xvec, h, thetavec);
    for(int i = 0; i < numtheta + 1; i++)
    {
      if(i == 0)
      {
        (gradient.slice(0)).col(curcol - 1) = k * gdmat.slice(0) * (qmattheta.slice(0)).col(curcol - 1);
      }
      else
      {
        (gradient.slice(i)).col(curcol - 1) = k * gdmat.slice(i) * (qmattheta.slice(0)).col(curcol - 1) + k * (gdmat.slice(0)) * (qmattheta.slice(i)).col(curcol - 1);
      }
    }
  }

return gradient;
}