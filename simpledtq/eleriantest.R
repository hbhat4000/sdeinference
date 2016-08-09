rm(list = ls(all = TRUE))

# load data
load('doublewell.RData')

# load necessary functions
library('sde')

drift <- function(t,x,theta)
{
  theta[1]*x*(theta[2] - x^2)
}

driftx <- function(t,x,theta)
{
  -3*x^2*theta[1] + theta[1]*theta[2]
}

driftxx <- function(t,x,theta)
{
  -6*x*theta[1]
}

diffusion <- function(t,x,theta)
{
  exp(theta[3])
}

diffusionx <- function(t,x,theta)
{
  return(0)
}

diffusionxx <- function(t,x,theta)
{
  return(0)
}

myh = 0.05
myk = myh^(0.75)
mybigm = ceiling(pi/(myk^1.5))
nsamples = 10000
thetamat = matrix(nrow=(nsamples+1),ncol=3)

source('mcmcstuff.R')

loglik = numeric(length=(nsamples+1))
logpost = numeric(length=(nsamples+1))
ar = numeric(length=nsamples)

ozakilik <- function(theta,h,k,bigm,littlet,data)
{
  lik = 0
  for (i in c(1:nrow(data)))
  {
    for (j in c(1:(ncol(data)-1)))
    {
      temp = dcElerian(x=data[i,j+1],t=littlet*j,x0=data[i,j],t0=littlet*(j-1),theta,d=drift,s=diffusion,sx=diffusionx,log=TRUE)
      lik = lik + temp
    }
  }
  return(lik)
}

loglik[1] = ozakilik(thetamat[1,],h=myh,k=myk,bigm=mybigm,littlet=1,data=xtraj)
# rawlik[rawlik <= 2.2e-16] = 0
# loglik[1] = sum(log(rawlik))
logpost[1] = loglik[1] + logprior(thetamat[1,])

# Metropolis algorithm
for (i in c(1:nsamples))
{
  thetastar = thetamat[i,] + propZ(n=1)
  thisloglik = ozakilik(thetastar,h=myh,k=myk,bigm=mybigm,littlet=1,data=xtraj)
  # rawlik[rawlik <= 2.2e-16] = 0
  # thisloglik = sum(log(rawlik))
  thislogpost = thisloglik + logprior(thetastar)
  ratio = exp(thislogpost - logpost[i])
  if (ratio > runif(n=1))
  {
    thetamat[i+1,] = thetastar
    loglik[i+1] = thisloglik
    logpost[i+1] = thislogpost
    ar[i] = 1
  }
  else
  {
    thetamat[i+1,] = thetamat[i,]
    loglik[i+1] = loglik[i]
    logpost[i+1] = logpost[i]
  }
  print(mean(ar[1:i]))
}

save.image(file="elerianresults.RData")


