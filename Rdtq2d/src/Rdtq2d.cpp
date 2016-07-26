#include "Rdtq2d.h"

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

// what should we pass in to solve the problem?
// 1) should pass in the full (x^r(t), y^r(t)) curve
//    together with the vector of times
// 2) nuvec should contain \nu_1, \nu_2, and \gamma(t)
// 3) \gamma(t) should be defined at the same vector of times as passed in;
//    could subsample so it's not so high-dimensional
// 4) should also pass in the local time at which function must be evaluated
// 5) then we interpolate and evaluate the whole thing

// assumptions:
// the matrix runner has three columns: time, x-coord, y-coord (of runner)

// f function for 2d
static inline vec fv(const vec& chase, const vec& gammavec, const mat& runner, double tval)
{
  // figure out time series time that is closest to t
  vec temptimes = tval - runner.col(0);

  unsigned int i;
  for (i=0; i<runner.n_rows; i++) if (temptimes(i) < 0) break;
  unsigned int ci = i-1;

  // interpolated runner's position
  vec ri = (runner((ci+1),span(1,2))*(tval - runner(ci,0)) + runner(ci,span(1,2))*(runner((ci+1),0) - tval)).t();
  double dt = runner((ci+1),0) - runner(ci,0);
  ri /= dt;

  // drift function as a vector
  vec drift = ri - chase;
  double denom = norm(drift);
  drift /= denom;
  drift *= gammavec(ci);
  return(drift);
}

// g function for 2d
static inline vec gv(const vec& chase, const vec& nuvec)
{
  return(nuvec);
}

mat dtq(const vec &nuvec, const vec &gammavec, const mat &runner, const mat &chaser, double h, int numsteps, double k, double yM)
{
  vec tvec = runner.col(0);
  vec C1 = chaser.col(1);
  vec C2 = chaser.col(2);

  double h12 = sqrt(h);
  int M = ceil(yM/k);
  int veclen = 2*M+1;
  int nr = veclen*veclen;
  int datapoints = C1.n_elem;
  vec xvec = k*linspace<vec>(-M,M,veclen);
  mat approxpdfvec(nr,datapoints);

/*
  vec chaser(2);
  chaser(0) = 0.123;
  chaser(1) = 0.246;
  vec ourcheck = fv(chaser, gammavec, runner, h);
*/

omp_set_num_threads(24);

#pragma omp parallel for
  for (int i = 0; i < datapoints; i++)
  {
    vec yvec(2);
    yvec(0) = C1(i);
    yvec(1) = C2(i);
    vec ftemp = fv(yvec, gammavec, runner, tvec(i));
    vec gtemp = gv(yvec, nuvec);
    vec v1 = gaussian_pdf_vec(xvec, C1(i) + ftemp(0)*h, h12*gtemp(0));
    vec v2 = gaussian_pdf_vec(xvec, C2(i) + ftemp(1)*h, h12*gtemp(1));

    approxpdfvec.col(i) = vectorise(kron(v1, v2.t()));
  }
#pragma omp barrier

  double supg = 1;
  int gamma = ceil(5*h12*supg/k);

  // loop over the timesteps
  for (int step=1; step<numsteps; step++)
  {
    mat approxpdfvecnew = zeros<mat>(nr,datapoints);
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
              vec yvec(2);
              yvec(0) = (ip-M)*k;
              yvec(1) = (jp-M)*k;
              vec sigma = h12*gv(yvec, nuvec);
              rowvec locbigg1(datapoints);
              rowvec locbigg2(datapoints);
              for (int tp=0; tp<datapoints; tp++)
              {
                double lefttime = runner(tp,0);
                vec mu = yvec + h*fv(yvec, gammavec, runner, (lefttime+step*h));
                locbigg1(tp) = gaussian_pdf(xvec(i), mu(0), sigma(0));
                locbigg2(tp) = gaussian_pdf(xvec(j), mu(1), sigma(1));
              }
              approxpdfvecnew.row(newspace) += k*k*(locbigg1 % locbigg2 % approxpdfvec.row(test));
            }
          }
        }
      }
    }
#pragma omp barrier
    approxpdfvec = approxpdfvecnew;
  }
  return approxpdfvec;
}

mat PDFcheck(const vec &nuvec, const vec &gammavec, const mat &runner, double h, double k, double yM)
{
  double h12 = sqrt(h);
  int M = ceil(yM/k);
  int veclen = 2*M+1;
  int nr = veclen*veclen;
  vec xvec = k*linspace<vec>(-M,M,veclen);

  double supg = 0.5;
  int gamma = ceil(5*h12*supg/k);
  int dimg = 2*gamma;

  mat Gnorm = zeros<mat>(veclen-dimg,veclen-dimg);

#pragma omp parallel for
  for (int newspace=0; newspace<nr; newspace++)
  {
    int j = floor(newspace/veclen); // subtract M to get math i, j \in [-M,M]
    int i = newspace - veclen*floor(newspace/veclen);      
    for (int ip=gamma; ip<(veclen-gamma); ip++)
    {
      for (int jp=gamma; jp<(veclen-gamma); jp++)
      {
        int test = jp*veclen + ip;
        if ((test >= 0) && (test < nr))
        {
          {
            vec yvec(2);
            yvec(0) = (ip-M)*k;
            yvec(1) = (jp-M)*k;
            vec mu = yvec + h*fv(yvec, gammavec, runner, h);
            vec sigma = h12*gv(yvec, nuvec);
            double locbigg1 = gaussian_pdf(xvec(i), mu(0), sigma(0));
            double locbigg2 = gaussian_pdf(xvec(j), mu(1), sigma(1));
            Gnorm(ip-gamma,jp-gamma) += k*k*locbigg1*locbigg2;
          }
        }
      }
    }
  }
#pragma omp barrier
  return Gnorm;
}

SEXP dtq2dCPP(SEXP s_nuvec, SEXP s_gammavec, SEXP s_runner, SEXP s_chaser, SEXP s_h, SEXP s_numsteps, SEXP s_k, SEXP s_yM)
{
    vec nuvec = Rcpp::as<arma::vec>(s_nuvec);
    vec gammavec = Rcpp::as<arma::vec>(s_gammavec);
    mat runner = Rcpp::as<arma::mat>(s_runner);
    mat chaser = Rcpp::as<arma::mat>(s_chaser);
    // cout << runner.n_cols << endl;
    // cout << runner.n_rows << endl;
    // cout << chaser.n_cols << endl;
    // cout << chaser.n_rows << endl;

    // cout << runner.col(0) << endl;
    // cout << runner.col(1) << endl;
    // cout << runner.col(2) << endl;
    // cout << chaser.col(0) << endl;
    // cout << chaser.col(1) << endl;
    // cout << chaser.col(2) << endl;

    double h = Rcpp::as<double>(s_h);
    int numsteps = Rcpp::as<int>(s_numsteps);
    double k = Rcpp::as<double>(s_k);
    double yM = Rcpp::as<double>(s_yM);

    mat mymat = dtq(nuvec, gammavec, runner, chaser, h, numsteps, k, yM);
    return Rcpp::wrap( mymat );
}


SEXP GCPP(SEXP s_nuvec, SEXP s_gammavec, SEXP s_runner, SEXP s_h, SEXP s_k, SEXP s_yM)
{
    vec nuvec = Rcpp::as<arma::vec>(s_nuvec);
    vec gammavec = Rcpp::as<arma::vec>(s_gammavec);
    mat runner = Rcpp::as<arma::mat>(s_runner);
    double h = Rcpp::as<double>(s_h);
    double k = Rcpp::as<double>(s_k);
    double yM = Rcpp::as<double>(s_yM);

    mat Gnorm = PDFcheck(nuvec, gammavec, runner, h, k, yM);
    return Rcpp::wrap( Gnorm );
}


