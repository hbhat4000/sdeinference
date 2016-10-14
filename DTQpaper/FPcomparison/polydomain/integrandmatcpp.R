require(Rcpp)
sourceCpp(code = '#include <Rcpp.h>
  #include <omp.h>
  using namespace Rcpp;
  // [[Rcpp::export]]
  NumericMatrix dnormpar2(NumericVector x, NumericVector mu, NumericVector sig)
  {
    double nc = 1/sqrt(2*PI);
    int n = x.size();
    int muSize = mu.size();
    int nn = n*muSize;
    NumericMatrix ret(n,muSize);
  omp_set_num_threads(12);
  #pragma omp parallel for
    for (int ii=0; ii<nn; ++ii)
    {
      int i = floor(ii/muSize);
      int j = ii % muSize;
      double s0 = sig(j);
      double x0 = x(i) - mu(j);
      ret(i,j) = exp(-x0*x0/(2*s0*s0))*nc/s0;
    }

    return ret;
  }'
)

integrandmat <- function(xvec,yvec,h,driftfun,difffun)
{
  drifty = driftfun(yvec)
  diffy = difffun(yvec)
  mymean = yvec+drifty*h
  mysd = diffy*sqrt(h)
  out = Matrix(dnormpar2(xvec,mymean,mysd))
  return(out)
}

