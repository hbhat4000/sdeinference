rm(list = ls(all = TRUE))

# load data
# load('doublewell.RData')

# load necessary functions
source('dtq.R')

myh = 0.05
myk = myh^(0.75)
mybigm = ceiling(pi/(myk^1.5))
nsamples = 10000
thetamat = matrix(nrow=(nsamples+1),ncol=3)

# source('mcmcstuff.R')
# define log prior
logprior <- function(z)
{
  return(sum(dnorm(x = z, mean = 0, sd = 100, log = TRUE)))
}
thetamat[1,] = c(1,1,1)

propZ <- function(n)
{
  return(rnorm(n, sd = 0.025))
}

loglik = numeric(length=(nsamples+1))
logpost = numeric(length=(nsamples+1))
ar = numeric(length=nsamples)

rawlik = dtq(thetamat[1,],h=myh,k=myk,bigm=mybigm,littlet=1,data=xtraj)
rawlik[rawlik <= 2.2e-16] = 0
loglik[1] = sum(log(rawlik))
logpost[1] = loglik[1] + logprior(thetamat[1,])

# Metropolis algorithm
for (i in c(1:nsamples))
{
  thetastar = thetamat[i,] + propZ(n=1)
  rawlik = dtq(thetastar,h=myh,k=myk,bigm=mybigm,littlet=1,data=xtraj)
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