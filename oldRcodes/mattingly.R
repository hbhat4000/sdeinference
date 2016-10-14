# clear memory
rm(list=ls(all=TRUE))

# drift and diffusion functions
driftfun <- function(x) { return(-x) }
difffun <- function(x) { return(rep(1,length(x))) }

integrandmat <- function(xvec,yvec,h,driftfun,difffun)
{
  X=replicate(length(yvec),xvec)
  Y=t(replicate(length(xvec),yvec))
  out = exp(-(X-Y-driftfun(Y)*h)^2/(2*difffun(Y)^2*h))/(difffun(Y)*sqrt(2*pi*h))
  return(out)
}

# simulation parameters
T = 1
s = 0.75
h = 0.025
init = 0
numsteps = ceiling(T/h)
k = h^s
yM = k*(pi/(k^2))
M = ceiling(yM/k)
xvec = k*c(-M:M)

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
theta = 0.5
a1 = alpha1(theta)
a2 = alpha2(theta)

# set up transition matrix, one column at a time
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
Aold = integrandmat(xvec,xvec,h,driftfun,difffun)

# pdf after one time step with Dirac \delta(x-init) initial condition
mymu = init + driftfun(init)*h
mysigma = abs(difffun(init))*sqrt(h)
phat = as.matrix(dnorm(x=xvec,mean=mymu,sd=mysigma))

# main iteration loop
for (i in c(2:numsteps)) phat = k*(A%*%phat)

# compare solutions
plot(xvec,phat)
truepdf = exp(-xvec^2/(1 - exp(-2*T)))/sqrt(pi*(1-exp(-2*T)))
lines(xvec,truepdf,col='red')

# L1 errs
# h = .025, 8.325528e-05
# h = .05, 0.0003507056
# h = .1, 0.001560497
# h = .25, 0.01284107
# h = .5, 0.07446556

hvec = c(0.025, 0.05, 0.1, 0.25, 0.5)
errvec = c(8.325528e-05, 0.0003507056, 0.001560497, 0.01284107, 0.07446556)
