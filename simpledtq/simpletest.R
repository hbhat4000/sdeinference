rm(list = ls(all = TRUE))

# load data
load('doublewell.RData')

# load necessary functions
source('dtq.R')

myh = 0.05
myk = myh^(0.75)
mybigm = ceiling(pi/(myk^1.5))
nsamples = 10000
thetamat = matrix(nrow=(nsamples+1),ncol=3)
thetamat[1,] = c(0.25, 4, 0)

propZ <- function(n)
{
  outmat = matrix(0,nrow=n,ncol=3)
  outmat[,1] = rnorm(n=n,mean=0,sd=0.02)
  outmat[,2] = rnorm(n=n,mean=0,sd=0.02)
  outmat[,3] = rnorm(n=n,mean=0,sd=0.01)
  return(outmat)
}

logprior <- function(theta)
{
  out = dnorm(theta[1],mean=0.5,sd=4,log=TRUE)
  out = out + dnorm(theta[2],mean=0.5,sd=4,log=TRUE)
  out = out + dnorm(theta[3],mean=0,sd=4,log=TRUE)
  return(out)
}

loglik = numeric(length=(nsamples+1))
logpost = numeric(length=(nsamples+1))
ar = numeric(length=nsamples)

rawlik = cdt(thetamat[1,],h=myh,k=myk,bigm=mybigm,littlet=1,data=xtraj)
rawlik[rawlik <= 2.2e-16] = 0
loglik[1] = sum(log(rawlik))
logpost[1] = loglik[1] + logprior(thetamat[1,])

# Metropolis algorithm
for (i in c(1:nsamples))
{
  thetastar = thetamat[i,] + propZ(n=1)
  rawlik = cdt(thetastar,h=myh,k=myk,bigm=mybigm,littlet=1,data=xtraj)
  rawlik[rawlik <= 2.2e-16] = 0
  thisloglik = sum(log(rawlik))
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



