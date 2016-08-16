# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data which has 2 lists, chaser and runner, each with t, x, y
chaserlist = list(NULL)
runnerlist = list(NULL)

# which deltat file do we use?
deltat = 0.4

xymax = 0
for (iii in c(1:1))
{
  fname = paste('fakedataswitch4/fakedata_dt_',deltat,'.RData',sep='')
  load(fname)
  chaserlist[[iii]] = matrix(unlist(chaser), ncol = 3)
  runnerlist[[iii]] = matrix(unlist(runner), ncol = 3)
  xymax = max(xymax,max(chaserlist[[iii]][,-1]))
  xymax = max(xymax,max(runnerlist[[iii]][,-1]))
}

# time increment from data
timeinc = (runnerlist[[1]])[2,1] - (runnerlist[[1]])[1,1]

# DTQ algorithm parameters
myh = timeinc/1
myk = timeinc/1 # 0.05*(myh)^(1.2)
xylimit = 1.5*xymax
gammalen = 8/myh

# number of DTQ steps
myns = floor(timeinc/myh)

# MCMC parameters
nsamples = 10000
burnin = 100
mcmcsteps  = nsamples + burnin
thetadim = 4
thetamat = matrix(0,nrow=(mcmcsteps+1),ncol=thetadim)
likelihoods = numeric(length=(mcmcsteps+1))
logposts = numeric(length=(mcmcsteps+1))
artrack = numeric(length=mcmcsteps)

# prior
meanvec = c(log(0.576222),log(0.576222),log(0.4),log(0.4))
logprior <- function(theta)
{
    out = dnorm(x=theta[1],mean=meanvec[1],sd=1,log=TRUE)
    out = out + dnorm(x=theta[2],mean=meanvec[2],sd=1,log=TRUE)
    out = out + dnorm(x=theta[3],mean=meanvec[3],sd=1,log=TRUE)
    out = out + dnorm(x=theta[4],mean=meanvec[4],sd=1,log=TRUE)
    return(out)
}

# proposal
proposalZ <- function(scaling=1)
{
    step = numeric(length=4)
    step[1] = rnorm(n=1,mean=0,sd=0.05/scaling)
    step[2] = rnorm(n=1,mean=0,sd=0.05/scaling)
    step[3] = rnorm(n=1,mean=0,sd=0.02/scaling)
    step[4] = rnorm(n=1,mean=0,sd=0.02/scaling)
    return(step)
}

# likelihood function
mylik <- function(theta)
{
    nuvec = exp(theta[3:4])
    gammavec = numeric(length=gammalen)
    gammavec[1:(gammalen/2)] = exp(theta[1])
    gammavec[(gammalen/2 + 1):gammalen] = exp(theta[2])
    den = Rdtq2d(nuvec, gammavec, runnerlist[[1]], chaserlist[[1]], myh, myns, myk, xylimit)
    loglik = sum(log(den[den >= 2.2e-16]))
    return(loglik)
}

# first guess
thetamat[1,] = meanvec
likelihoods[1] = mylik(thetamat[1,])
logposts[1] = likelihoods[1] + logprior(thetamat[1,])

# main MCMC loop
for (i in c(1:mcmcsteps))
{
  thetastar = thetamat[i,] + proposalZ(scaling=0.15)
  likelihood = mylik(thetastar)
  logpost = likelihood + logprior(thetastar)

  # accept/reject step
  rho = runif(n=1)
  ratio = exp(logpost - logposts[i])
  if (ratio > rho)
  {
    # accept
    artrack[i] = 1
    thetamat[(i+1),] = thetastar
    likelihoods[i+1] = likelihood
    logposts[i+1] = logpost
    print(paste("Accepted",paste(exp(thetastar),collapse=', ')))
  }
  else
  {
    # reject
    artrack[i] = 0
    thetamat[(i+1),] = thetamat[i,]
    likelihoods[i+1] = likelihoods[i]
    logposts[i+1] = logposts[i]
  }
  print(c(i,mean(artrack[1:i])))
}

thetasamp = thetamat[-c(1:burnin),]

fname = paste("mcmcout_",deltat,"_1.RData",sep='')
save.image(file=fname)



