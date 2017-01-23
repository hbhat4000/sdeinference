#include <iostream>
#include <armadillo>
#include <cmath>

// Gaussian functions
double gausspdfscalar(double x, double mu, double sigma)
{
  double u = (x - mu) / fabs(sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u * u) / 2);
  return p; 
}

arma::vec gausspdf(const arma::vec &x, double mu, double sigma)
{
  arma::vec u = (x - mu) / fabs(sigma);
  arma::vec p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-(u % u) / 2);
  return p; 
}

arma::vec gausspdf(double x, const arma::vec &mu, const arma::vec &sigma)
{
  arma::vec u = (x - mu) / abs(sigma);
  arma::vec p = (1 / (sqrt (2 * M_PI) * abs (sigma))) % exp (-(u % u) / 2);
  return p; 
}

// the following typedef is intended to be used for user-defined
// drift and diffusion functions; note that the function depends both
// on a spatial location x and the parameter vector theta
typedef double (*funcPtr)(const double& x, const arma::vec& theta);

// same as above except that the function returns a vector; this is
// suitable for gradients of f and g
typedef arma::vec (*gradPtr)(const double& x, const arma::vec& theta);

class dtq {
  // private data
  int spi; // DTQ steps per interval
  double myh; // DTQ time step
  double myk; // DTQ grid spacing
  int bigm; // DTQ grid goes from -bigm*myk to bigm*myk
  funcPtr f; // DTQ drift function 
  funcPtr g; // DTQ diffusion function

  double deltat; // time step of the data
  arma::vec curtheta; // current value of theta

  gradPtr gradf; // gradient of DTQ drift function w.r.t. theta
  gradPtr gradg; // gradient of DTQ diffusion function w.r.t. theta

  // flags so that we know whether something exists yet
  bool haveMyh = false;
  bool haveProp = false;
  bool haveData = false;
  bool haveGradfg = false;
  bool haveLoglik = false;
  bool haveGradloglik = false;

  // for next two, we use pointers to save memory
  arma::vec *tvec; // vector of times at which data was taken
                   // think of this as a column vector
  arma::mat *odata; // observations at the times given by tvec
                    // think of this as several columns, one per time series
                    // odata->n_rows should equal tvec->n_elem
  int ltvec;
  int numts;

  arma::sp_mat prop; // the propagator matrix (sparse)
  int ylen;
  arma::vec yvec;
  arma::vec fy;
  arma::vec gy;

  arma::vec loglikvec; // log likelihood at each consecutive pair of times
  double loglik; // total log likelihood
  arma::vec gradloglik; // gradient of total log likelihood w.r.t. theta

  public:
    // constructors
    // minimal constructor for likelihood calculation
    dtq(int inspi, double ink, int inm, funcPtr inf, funcPtr ing, arma::vec &theta) : 
      spi(inspi), myk(ink), bigm(inm), f(inf), g(ing), curtheta(theta) {}

    // minimal constructor needed to solve one forward problem
    dtq(double inmyh, double ink, int inm, funcPtr inf, funcPtr ing, arma::vec &theta) : 
      myh(inmyh), myk(ink), bigm(inm), f(inf), g(ing), curtheta(theta) { haveMyh = true; }

    // retrieve mandatory bits of data
    double getH(void) { return myh; }
    double getK(void) { return myk; }
    int getBigm(void) { return bigm; }

    // retrieve other bits of data
    int getLtvec(void);
    int getNumts(void);

    // setting gradients
    int setGrads(gradPtr, gradPtr);

    // compute the propagator
    int compProp(void);

    // setting/loading data
    int setData(arma::vec *, arma::mat *);

    // compute log likelihood and gradient
    int compLL(void);
    int compGrad(void);

    // retrieve log likelihood and gradient
    double getLL();
    int getGrad(arma::vec &);
};

int dtq::setGrads(gradPtr ingradf, gradPtr ingradg)
{
  gradf = ingradf;
  gradg = ingradg; 
  haveGradfg = true;
  return 0;
}

int dtq::compProp(void)
{
  // need myh to proceed
  if (! haveMyh) return 1;

  // fundamental size
  ylen = 2*bigm+1;

  // start computing the propagator!
  prop = arma::sp_mat(ylen,ylen);

  // create yvec
  yvec = arma::zeros(ylen);
  for (int i=-bigm; i<=bigm; i++)
    yvec(bigm+i) = i*myk;

  // apply f and g to yvec
  fy = arma::zeros(ylen);
  gy = arma::zeros(ylen);
  for (int i=-bigm; i<=bigm; i++)
  {
    fy(bigm+i) = (*f)(yvec(bigm+i),curtheta);
    gy(bigm+i) = (*g)(yvec(bigm+i),curtheta);
  }

  // normalization "constant"
  arma::vec c0mod = 1.0/(sqrt(2.0*(arma::datum::pi)*myh)*gy);

  // variance
  arma::vec myvar = gy*gy*myh;

  // compute and set main diagonal
  arma::vec maindiag = exp(-(myh/2.0)*(fy%fy)/(gy*gy)) % c0mod;
  prop.diag() = maindiag;

  // superdiagonals
  bool done = false;
  int curdiag = 1;
  double mytol = 2.0e-16;
  double refsum = arma::sum(arma::abs(maindiag))*myk;
  while (! done)
  {
    arma::vec mymean = curdiag*myk + fy*myh;
    arma::vec thisdiag = exp(-mymean%mymean/(2.0*myvar))%c0mod;
    thisdiag = thisdiag.tail(ylen - curdiag);
    double thissum = arma::sum(arma::abs(thisdiag));
    if ((curdiag == 1) || (thissum > mytol*refsum))
    {
      prop.diag(curdiag) = thisdiag;
    }
    else done = true;
  }
  int maxdiag = curdiag;
  for (curdiag=1; curdiag<=maxdiag; curdiag++)
  {
    arma::vec mymean = -curdiag*myk + fy*myh;
    arma::vec thisdiag = exp(-mymean%mymean/(2.0*myvar))%c0mod;
    thisdiag = thisdiag.head(ylen - curdiag);
    prop.diag(-curdiag) = thisdiag;
  }
  
  haveProp = true;
  return 0;
}

int dtq::setData(arma::vec *intvec, arma::mat *inodata)
{
  int lintvec = (int) intvec->n_elem;
  int inodatac = (int) odata->n_cols;
  int inodatar = (int) odata->n_rows;
  // check whether we have at least two points in time,
  // whether number of rows of data equals number of points in time,
  // and whether we have at least one sample path
  if ((lintvec >= 2) && (lintvec==inodatar) && (inodatac>=1))
  {
    tvec = intvec;
    odata = inodata;
    deltat = (*tvec)(1) - (*tvec)(0); // assume equispaced data
    if (haveMyh)
    {
      spi = ceil(deltat/myh);      // or set # of DTQ steps per interval
    }
    else
    {
      myh = deltat/spi;   // set DTQ time step
      haveMyh = true;
    }
    ltvec = lintvec;
    numts = inodatac;
    haveData = true;
    return 0;
  }
  else
  {
    tvec = NULL;
    odata = NULL;
    haveData = false;
    return 1;
  }
}

int dtq::getLtvec(void)
{
  if (haveData)
    return ltvec;
  else
    return 1;
}

int dtq::getNumts(void)
{
  if (haveData)
    return numts;
  else
    return 1;
}


int dtq::compLL(void)
{
  // remember, everything here is for equispaced data
  // we'll save the non-equispaced case for our scala + spark code :)
  if ((! haveData) || (! haveMyh)) return 1;
  if (spi<1) return 1;

  loglikvec = arma::zeros(ltvec-1);

  if (spi==1) // special case
  {
    for (int i=0; i<(ltvec-1); i++)
    {
      for (int j=0; j<numts; j++)
      {
        double xi = (*odata)(i,j);
        double xip1 = (*odata)(i+1,j);
        double mu = xi + ((*f)(xi,curtheta))*myh;
        double myh12 = sqrt(myh);
        double sig = ((*g)(xi,curtheta))*myh12;
        loglikvec(i) += gausspdfscalar(xip1,mu,sig);
      }
    }
  }
  else
  {
    // build the big matrix of initial conditions
    arma::mat dtqmat = arma::zeros(ylen,(ltvec-1));
    double myh12 = sqrt(myh);
    for (int i=0; i<(ltvec-1); i++)
    {
      // go through each particular initial condition at this time
      // and make a Gaussian
      for (int j=0; j<numts; j++)
      {
        double xi = (*odata)(i,j);
        double mu = xi + ((*f)(xi,curtheta))*myh;
        double sig = ((*g)(xi,curtheta))*myh12;
        dtqmat.col(i) += gausspdf(yvec,mu,sig);
      }
      dtqmat.col(i) = dtqmat.col(i) / numts;
    }

    // propagate this forward in time by (spi-2) steps
    if (spi >= 3)
      for (int i=1; i<=(spi-2); i++)
        dtqmat = myk * prop * dtqmat;

    // now multiply on the left by the Gamma vectors
    arma::vec muvec = yvec + fy*myh;
    arma::vec sigvec = gy*myh12;
    for (int i=0; i<(ltvec-1); i++)
    {
      for (int j=0; j<numts; j++)
      {
        arma::vec gammavec = myk*gausspdf((*odata)(i+1,j),muvec,sigvec);
        loglikvec(i) += arma::dot(gammavec,dtqmat.col(i));
      }
    }
  }
  haveLoglik = true;
  return 0;
}

int main(void)
{
  // create fake data and compute log likelihood

  arma::vec test = arma::zeros(3);
  arma::vec *testptr;
  testptr = &test;
  test(0) = 0.5;
  test(1) = 2.0;
  test(2) = -1.0;
  std::cout << (*testptr)(1) << '\n';
}


