# anderson & mattingly functions and parameters
alpha1 <- function(th)
{
  return(1/(2*th*(1-th)))
}
alpha2 <- function(th)
{
  num = (1-th)^2 + th^2
  denom = 2*th*(1-th)
  return(num/denom)
}
rho <- function(x)
{
  v = x
  v[x<0] = 0
  return(v)
}

integrandmat <- function(xvec,yvec,h,driftfun,difffun)
{
  theta = 0.5
  a1 = alpha1(theta)
  a2 = alpha2(theta)

  A = matrix(0,nrow=(2*M+1),ncol=(2*M+1))
  for (i in c(1:(2*M+1)))
  {
    xjm1 = xvec[i]
    mu1 = xjm1 + driftfun(xjm1)*theta*h
    sig1 = abs(difffun(xjm1))*sqrt(theta*h)
    pvec = dnorm(x=xvec,mean=mu1,sd=sig1)
    xjmat = replicate((2*M+1),xvec)
    xstarmat = t(xjmat)
    mu2 = xstarmat + (a1*driftfun(xstarmat) - a2*driftfun(xjm1))*(1-theta)*h
    sig2 = sqrt(rho(a1*difffun(xstarmat)^2 - a2*difffun(xjm1)^2))*sqrt((1-theta)*h)
    pmat = exp(-(xjmat-mu2)^2/(2*sig2*sig2))/(sig2*sqrt(2*pi))
    A[,i] = k*(pmat %*% pvec)
  }
  return(A)
}


