#include "RgdtqSUM.h"

#include <RcppArmadillo.h>
#include <cmath>

using namespace std;
using namespace arma;

/* Linear case */
// drift function
static inline double f(const double y, const vec& thetavec)
{
  return(thetavec(0) * (thetavec(1) - y));
}

static inline double g(const double y, const vec& thetavec)
{
  return(thetavec(2));
}

// Gaussian PDF and Gradients
// takes as input 2 points, x and y
// returns a vector of computed gaussian pdf and derivatives wrt thetavec
static inline vec gaussian_pdf(const double x, const double y, const double h, const vec& thetavec)
{
  double driftval = f(y, thetavec);
  double sigma = abs(g(y, thetavec));
  double sqsigma = pow(sigma, 2);
  double part1 = x - y - driftval*h;
  vec gaussian = zeros<vec>(thetavec.n_elem + 1);

  double G = exp (- (part1 * part1) / (2 * h * sqsigma)) / (sqrt (2 * M_PI * h) * sigma);
  double dGdf = G * part1 / sqsigma;
  double dGdg = G * (-1 / sigma + (part1 * part1)/(sigma * sqsigma * h));

  gaussian[0] = G;
  gaussian[1] = dGdf * (driftval / thetavec(0));
  gaussian[2] = dGdf * thetavec(0);
  gaussian[3] = dGdg;

  return gaussian;
}

cube gDTQ(const vec& thetavec, const double h, const double k, const int M, const int littlet, const mat &init_data)
{
  int numsteps = ceil(littlet/h);

  int veclen = 2*M + 1;
  int datapoints = init_data.n_cols;
  vec xvec = k*(linspace<vec>(-M, M, veclen));
  int numtheta = thetavec.n_elem;
  // double supg = 0.5;
  // int gamma = ceil(5 * sqrt(h) * supg / k);
  // int dimg = 2 * gamma;

/**** NEW SEGMENT ****/
  cube DSUM = zeros<cube>(veclen, veclen, numtheta + 1);
  cube qmatthetaSUM = zeros<cube>(veclen, datapoints - 1, numtheta + 1);

  for(int i = 0; i < veclen; i++)
  {
    for(int j = 0; j < veclen; j++)
    {
      int x = xvec(i);
      int y = xvec(j);

      // gaussian_pdf returns a vector of length numtheta + 1
      // the subcube is selected as Dsum(i, j, :)
      DSUM.subcube(span(i), span(j), span()) = gaussian_pdf(x, y, h, thetavec);

      // // Alternate way which is equivalent 
      // vec Dvec = gaussian_pdf(x, y, h, thetavec);
      // DSUM.slice(0)(i, j) = Dvec(0);
      // DSUM.slice(1)(i, j) = Dvec(1);
      // DSUM.slice(2)(i, j) = Dvec(2);
      // DSUM.slice(3)(i, j) = Dvec(3);
    }
  }
  cout << "Dimensions of cube DSUM: " << DSUM.n_rows << ", " << DSUM.n_cols << ", " << DSUM.n_slices << endl;
  cout << "Sum of cube DSUM: " << accu(DSUM) << endl;
  cout << "3*3 matrix for objective function: " << DSUM.subcube(0, 0, 0, 2, 2, 3) << endl;

/**** NEW SEGMENT ****/  
  cube initderivsSUM = zeros<cube>(veclen, init_data.n_rows, numtheta + 1);
  for(int curcol = 0; curcol < datapoints - 1; curcol++)
  {
    for(int i = 0; i < veclen; i++)
    {
      for(int j = 0; j < datapoints - 1; j++)
      {
        int x = xvec(i);
        int y = init_data(j, curcol);

        // gaussian_pdf returns a vector of length numtheta + 1
        // the subcube is selected as Dsum(i, j, :)
        initderivsSUM.subcube(span(i), span(j), span()) = gaussian_pdf(x, y, h, thetavec);
      }
    }
  }
  cout << "Dimensions of cube initderivsSUM:" << initderivsSUM.n_rows << ", " << initderivsSUM.n_cols << ", " << initderivsSUM.n_slices << endl;
  cout << "Sum of cube initderivsSUM: " << accu(initderivsSUM) << endl;

  for(int curcol = 0; curcol < datapoints - 1; curcol++)
  {
    for(int thetavals = 0; thetavals < numtheta + 1; thetavals++)
    {
      (qmatthetaSUM.slice(thetavals)).col(curcol) = mean(initderivsSUM.slice(thetavals), 1);
    }
  }
  cout << "Dimensions of cube qmatthetaSUM:" << qmatthetaSUM.n_rows << ", " << qmatthetaSUM.n_cols << ", " << qmatthetaSUM.n_slices << endl;
  cout << "Sum of cube qmatthetaSUM: " << accu(qmatthetaSUM) << endl;

/**** OLD SEGMENT reused because gaussian_pdf function not called****/
  int startstep = 1;
  for(int curstep = startstep; curstep < numsteps - 1; curstep++)
  {
    qmatthetaSUM.slice(0) = k * DSUM.slice(0) * qmatthetaSUM.slice(0);

    for(int i = 1; i < numtheta + 1; i++)
    {
      qmatthetaSUM.slice(i) = k * DSUM.slice(0) * qmatthetaSUM.slice(i) + k * DSUM.slice(i) * qmatthetaSUM.slice(0);
    }
  }

/**** NEW SEGMENT ****/
  cube gdmatSUM = zeros<cube>(init_data.n_rows, veclen, numtheta + 1);
  for(int curcol = 1; curcol < datapoints; curcol++)
  {
    for(int i = 0; i < datapoints - 1; i++)
    {
      for(int j = 0; j < veclen; j++)
      {
        int x = init_data(i, curcol);
        int y = xvec(j);

        // gaussian_pdf returns a vector of length numtheta + 1
        // the subcube is selected as Dsum(i, j, :)
        gdmatSUM.subcube(span(i), span(j), span()) = gaussian_pdf(x, y, h, thetavec);
      }
    }
  }
  cout << "Dimensions of cube gdmatSUM:" << gdmatSUM.n_rows << ", " << gdmatSUM.n_cols << ", " << gdmatSUM.n_slices << endl;
  cout << "Sum of cube gdmatSUM: " << accu(gdmatSUM) << endl;

/**** OLD SEGMENT reused after removing the gaussian_pdf function call****/  
// gradient.slice(0) = likelihood
// gdmat.slice(0) = gammamat
  cube gradientSUM = zeros<cube>(init_data.n_rows, datapoints - 1, numtheta + 1);
  for(int curcol = 1; curcol < datapoints; curcol++)
  {
    // cube gdmat = Dtheta(init_data.col(curcol), xvec, h, thetavec);
    (gradientSUM.slice(0)).col(curcol - 1) = k * gdmatSUM.slice(0) * (qmatthetaSUM.slice(0)).col(curcol - 1);

    for(int i = 1; i < numtheta + 1; i++)
    {
      (gradientSUM.slice(i)).col(curcol - 1) = k * gdmatSUM.slice(i) * (qmatthetaSUM.slice(0)).col(curcol - 1) + k * (gdmatSUM.slice(0)) * (qmatthetaSUM.slice(i)).col(curcol - 1);
    }
  }
  cout << "Dimensions of cube gradientSUM:" << gradientSUM.n_rows << ", " << gradientSUM.n_cols << ", " << gradientSUM.n_slices << endl;

  // // to check that the right values are getting passed to R
  cout << "objective: " << -accu(log(gradientSUM.slice(0))) << endl;
  cout << "gradientSUM 1: " << -accu(gradientSUM.slice(1)/gradientSUM.slice(0)) << endl;
  cout << "gradientSUM 2: " << -accu(gradientSUM.slice(2)/gradientSUM.slice(0)) << endl;
  cout << "gradientSUM 3: " << -accu(gradientSUM.slice(3)/gradientSUM.slice(0)) << endl;

return gradientSUM;
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
