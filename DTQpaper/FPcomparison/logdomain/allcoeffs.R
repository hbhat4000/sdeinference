### Example 1:
### dX_t = 1/2 X_t dt + sqrt(1 + X_t^2)dW_t

examples <- list(NULL)

examples[[1]]=list(NULL)
examples[[1]]$drift <- function(y) 
{
  return(y/2)
}

examples[[1]]$diff <- function(y)
{
  return(sqrt(1+y^2))
}

examples[[1]]$exact <- function(xvec,T)
{
  exact = exp(-(asinh(xvec)^2/(2*T)))/(sqrt(2*pi*T*(xvec^2+1)))
  return(exact)
}

### Example 2:
### dX_t = -1/2 tanh(X_t)sech^2(X_t)dt+sech(X_t)dW_t

examples[[2]]=list(NULL)
examples[[2]]$drift <- function(y)
{
  return(-tanh(y)*(1/cosh(y)^2)/2)
}

examples[[2]]$diff <- function(y)
{
  return(1/cosh(y))
}

examples[[2]]$exact <- function(xvec,T)
{
  exact = cosh(xvec)*exp(-(sinh(xvec)^2/(2*T)))/sqrt(2*pi*T)
  return(exact)
}

### Example 3:
### dX_t = -X_t dt + dW_t
examples[[3]]=list(NULL)
examples[[3]]$drift <- function(y)
{
  return(-y)
}

examples[[3]]$diff <- function(y)
{
  return(rep(1,length(y)))
}

examples[[3]]$exact <- function(xvec,T)
{
  normcoef = (pi*(1 - exp(-2*T)))^(-0.5)
  return(normcoef*exp(-xvec^2/(1 - exp(-2*T))))
}

### Example 4:
### dY_t = (-sqrt(1+Y_t^2) sinh^(-1)(Y_t) + (1/2) Y_t) dt + sqrt(1+Y_t^2) dW_t
examples[[4]]=list(NULL)
examples[[4]]$drift <- function(y)
{
  return(-sqrt(1+y^2)*asinh(y) + 0.5*y)
}

examples[[4]]$diff <- function(y)
{
  return(sqrt(1+y^2))
}

examples[[4]]$exact <- function(y,T)
{
  normcoef = (pi*(1 - exp(-2*T))*(1 + y^2))^(-0.5)
  return(normcoef*exp(-(asinh(y))^2/(1 - exp(-2*T))))
}

### Example 5:
### dX_t = (1/2) X_t dt + \sqrt{1 + X_t^2} dW_t
examples[[5]]=list(NULL)
examples[[5]]$drift <- function(y)
{
  return(0.5*y)
}
examples[[5]]$diff <- function(y)
{
  return(sqrt(1 + y^2))
}
examples[[5]]$sampler <- function(y)
{
  # n = number of samples
  # ft = final time
  # x0 = initial condition
  w = rnorm(n=y$n,mean=0,sd=sqrt(y$ft))
  return(sinh(w + asinh(y$x0)))
}
examples[[5]]$exact <- function(yin,T)
{
  fac1 = (1 + yin^2)^(-1/2)
  fac2 = dnorm(x=asinh(yin),mean=0,sd=sqrt(T))

#  # note: T is assumed to be 1
#  load('./sw_x_5.RData')
#  load('./sw_bigden_5.RData')
#  z = approx(x=mydenx,y=bigden,xout=yin)$y
#  k = yin[2]-yin[1]
#  zsum = k*sum(z)
#  z = z/zsum
#  print(k*sum(abs(z-fac1*fac2)))

  return(fac1*fac2)
}

### Example 6:
### dX_t = 
examples[[6]]=list(NULL)
examples[[6]]$drift <- function(y)
{
  return(-sin(y)*cos(y)^3)
}
examples[[6]]$diff <- function(y)
{
  return(cos(y)^2)
}
examples[[6]]$sampler <- function(y)
{
  w = rnorm(n=y$n,mean=0,sd=sqrt(y$ft))
  return(atan(w + tan(y$x0)))
}
#examples[[6]]$exact <- function(yin,T)
#{
#  # note: T is assumed to be 1
#  load('./sw_x_6.RData')
#  load('./sw_bigden_6.RData')
#  z = approx(x=mydenx,y=bigden,xout=yin)$y
#  zsum = (yin[2]-yin[1])*sum(z)
#  z = z/zsum
#  return(z)
#}

### Example 7:
### dX_t = 
examples[[7]]=list(NULL)
examples[[7]]$drift <- function(y)
{
  return(-0.5*tanh(y)*sech(y)^2)
}
examples[[7]]$diff <- function(y)
{
  return(sech(y))
}
examples[[7]]$sampler <- function(y)
{
  w = rnorm(n=y$n,mean=0,sd=sqrt(y$ft))
  return(asinh(w + sinh(y$x0)))
}
#examples[[7]]$exact <- function(yin,T)
#{
## note: T is assumed to be 1
#  load('./sw_x_7.RData')
#  load('./sw_bigden_7.RData')
#  z = approx(x=mydenx,y=bigden,xout=yin)$y
#  zsum = (yin[2]-yin[1])*sum(z)
#  z = z/zsum
#  return(z)
#}

### Example 8:
### dX_t = 
examples[[8]]=list(NULL)
examples[[8]]$drift <- function(y)
{
  return(0.5*y + sqrt(y^2 + 1))
}
examples[[8]]$diff <- function(y)
{
  return(sqrt(y^2+1))
}
examples[[8]]$sampler <- function(y)
{
  w = rnorm(n=y$n,mean=0,sd=sqrt(y$ft))
  return(sinh(y$ft + w + asinh(y$x0)))
}
examples[[8]]$exact <- function(xvec,T)
{
  # note: T is assumed to be 1
  # load('./sw_x_8.RData')
  # load('./sw_bigden_8.RData')
  # z = approx(x=mydenx,y=bigden,xout=xvec)$y
  # zsum = (yin[2]-yin[1])*sum(z)
  # z = z/zsum
  # return(z)
  transx = asinh(xvec) - T
  prefac = (1 + xvec^2)^(-1/2)
  z = prefac*dnorm(x=transx)
  return(z)
}


### Example 9:
### dX_t = -tanh(X_t) dt + dW_t
examples[[9]]=list(NULL)
examples[[9]]$drift <- function(y)
{
  return(-tanh(y))
}

examples[[9]]$diff <- function(y)
{
  return(rep(1,length(y)))
}

examples[[9]]$exact <- function(xvec,T)
{
#  load('./sw_x_9.RData')
#  load('./sw_bigden_9.RData')
#  z = approx(x=mydenx,y=bigden,xout=xvec)$y
#  zsum = (xvec[2]-xvec[1])*sum(z)
#  z = z/zsum
#  return(z)
  # empdf = density(emsol,from=-min(xvec),to=max(xvec),n=length(xvec),bw='nrd')
  # return(empdf$y)
}

examples[[9]]$sampler <- function(y)
{
  x = rep(y$x0,times=y$n)
  nsteps = ceiling(y$ft/y$h)
  for (i in c(1:nsteps))
    x = x - tanh(x)*y$h + sqrt(y$h)*rnorm(mean=0,sd=1,n=y$n)

  return(x)
}

