examples <- list(NULL)

### Example 1:
### dX_t = -X_t dt + dW_t
examples[[1]]=list(NULL)
examples[[1]]$drift <- function(y)
{
  return(-y)
}
examples[[1]]$diff <- function(y)
{
  return(rep(1,length(y)))
}
examples[[1]]$exact <- function(xvec,T)
{
  normcoef = (pi*(1 - exp(-2*T)))^(-0.5)
  return(normcoef*exp(-xvec^2/(1 - exp(-2*T))))
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
### dX_t = -(sin(X_t) cos^3(X_t)) dt + (cos^2(X_t)) dW_t
examples[[3]]=list(NULL)
examples[[3]]$drift <- function(y)
{
  return(-sin(y)*cos(y)^3)
}
examples[[3]]$diff <- function(y)
{
  return(cos(y)^2)
}
examples[[3]]$sampler <- function(y)
{
  w = rnorm(n=y$n,mean=0,sd=sqrt(y$ft))
  return(atan(w + tan(y$x0)))
}
examples[[3]]$exact <- function(yin,T)
{
  z = exp(-tan(yin)^2/(2*T))
  z = z/(cos(yin)^2)
  z = z/sqrt(2*pi*T)
  return(z)
  # note: T is assumed to be 1
  # load('./sw_x_3.RData')
  # load('./sw_bigden_3.RData')
  # z = approx(x=mydenx,y=bigden,xout=yin)$y
  # zsum = (yin[2]-yin[1])*sum(z)
  # z = z/zsum
  # return(z)
}

### Example 4:
### dX_t = ((1/2) X_t + sqrt(1+X_t^2)) dt + sqrt(1+X_t^2) dW_t
examples[[4]]=list(NULL)
examples[[4]]$drift <- function(y)
{
  return(0.5*y + sqrt(y^2 + 1))
}
examples[[4]]$diff <- function(y)
{
  return(sqrt(y^2+1))
}
examples[[4]]$sampler <- function(y)
{
  w = rnorm(n=y$n,mean=0,sd=sqrt(y$ft))
  return(sinh(y$ft + w + asinh(y$x0)))
}
examples[[4]]$exact <- function(xvec,T)
{
  transx = asinh(xvec) - T
  prefac = (1 + xvec^2)^(-1/2)
  z = prefac*dnorm(x=transx)
  return(z)
}

### Example 5:
### dX_t = 1/2 X_t dt + sqrt(1 + X_t^2)dW_t
examples[[5]]=list(NULL)
examples[[5]]$drift <- function(y) 
{
  return(y/2)
}
examples[[5]]$diff <- function(y)
{
  return(sqrt(1+y^2))
}
examples[[5]]$exact <- function(xvec,T)
{
  exact = exp(-(asinh(xvec)^2/(2*T)))/(sqrt(2*pi*T*(xvec^2+1)))
  return(exact)
}

### Example 6:
### dY_t = (-sqrt(1+Y_t^2) sinh^(-1)(Y_t) + (1/2) Y_t) dt + sqrt(1+Y_t^2) dW_t
examples[[6]]=list(NULL)
examples[[6]]$drift <- function(y)
{
  return(-sqrt(1+y^2)*asinh(y) + 0.5*y)
}
examples[[6]]$diff <- function(y)
{
  return(sqrt(1+y^2))
}
examples[[6]]$exact <- function(y,T)
{
  normcoef = (pi*(1 - exp(-2*T))*(1 + y^2))^(-0.5)
  return(normcoef*exp(-(asinh(y))^2/(1 - exp(-2*T))))
}

