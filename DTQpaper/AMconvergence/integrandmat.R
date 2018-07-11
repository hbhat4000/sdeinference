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

integrandmat <- function(xvec,k,h,driftfun,difffun)
{
  theta = 0.5
  a1 = alpha1(theta)
  a2 = alpha2(theta)

  twoMplus1 = length(xvec)
  A = matrix(0,nrow=twoMplus1,ncol=twoMplus1)
  xjmat = replicate(twoMplus1,xvec)
  xstarmat = t(xjmat)

  driftxvec = driftfun(xvec)
  diffxvec = difffun(xvec)

  drxstar = t(replicate(twoMplus1,driftxvec))
  dixstar2 = (t(replicate(twoMplus1,diffxvec)))^2
  
  # drxstar = driftfun(xstarmat)
  # dixstar = difffun(xstarmat)

  for (i in c(1:twoMplus1))
  {
    xjm1 = xvec[i]
    mu1 = xjm1 + driftxvec[i]*theta*h
    sig1 = abs(diffxvec[i])*sqrt(theta*h)
    pvec = dnorm(x=xvec,mean=mu1,sd=sig1)

    mu2 = xstarmat + (a1*drxstar - a2*driftxvec[i])*(1-theta)*h
    sig2 = sqrt(rho(a1*dixstar2 - a2*diffxvec[i]^2))*sqrt((1-theta)*h)
    pmat = dnorm(x=xjmat,mean=mu2,sd=sig2)
    pmat[abs(pmat)==Inf] = 0

    # pmat = exp(-(xjmat-mu2)^2/(2*sig2*sig2))/(sig2*sqrt(2*pi))
    temp = k*(pmat %*% pvec)
    temp[is.na(temp)] = 0
    A[,i] = temp
  }
  return(A)
}


