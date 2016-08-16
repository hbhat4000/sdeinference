# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load nba data
load('./nbaspatial/traj.RData')
traj = newdf

# rescale so that everything is between -0.5 and 0.5
traj[,2] = traj[,2]/94
traj[,3] = traj[,3]/50
traj[,4] = traj[,4]/94
traj[,5] = traj[,5]/50
traj[,2:5] = traj[,2:5] - 0.5

# make runner and chaser matrices
runnermat = as.matrix(traj[,1:3])
chasermat = as.matrix(cbind(traj[,1],traj[,4:5]))

# time increment from data
timeinc = traj[2,1] - traj[1,1]

# DTQ algorithm parameters
ndiv = 4
myh = timeinc/ndiv
myk = 2*myh
xylimit = 1.5*0.5
gammalen = 32
gammaveclen = traj[nrow(traj),1]/myh

# number of DTQ steps
myns = floor(timeinc/myh)

# MCMC parameters
nsamples = 1000
burnin = 100
mcmcsteps  = nsamples + burnin
thetadim = gammalen + 2
thetamat = matrix(0,nrow=(mcmcsteps+1),ncol=thetadim)
likelihoods = numeric(length=(mcmcsteps+1))
logposts = numeric(length=(mcmcsteps+1))
artrack = numeric(length=mcmcsteps)

# prior
# for gamma prior let's use the average speed of the chaser
x1 = diff(traj[,4])
x2 = diff(traj[,5])
meanspeed = mean(sqrt(x1^2 + x2^2)/timeinc)
meanvec = c(rep(log(meanspeed),gammalen),rep(log(0.4),2))
logprior <- function(theta)
{
    out = dnorm(x=theta,mean=meanvec,sd=1,log=TRUE)
    return(sum(out))
}

# proposal
proposalZ <- function(scaling=1)
{
    step = numeric(length=thetadim)
    step[1:gammalen] = rnorm(n=gammalen,mean=0,sd=0.05/scaling)
    step[(gammalen+1):(gammalen+2)] = rnorm(n=2,mean=0,sd=0.02/scaling)
    return(step)
}

# likelihood function
mylik <- function(theta)
{
    nuvec = exp(theta[(gammalen+1):(gammalen+2)])
    gammavec = numeric(length=gammaveclen)
    chunklen = floor(gammaveclen/gammalen)
    for (j in c(1:(gammalen-1)))
        gammavec[((j-1)*chunklen+1):(j*chunklen)] = exp(theta[j])

    gammavec[((gammalen-1)*chunklen+1):gammalen] = exp(theta[gammalen])

    den = Rdtq2d(nuvec, gammavec, runnermat, chasermat, myh, myns, myk, xylimit)
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
  thetastar = thetamat[i,] + proposalZ(scaling=0.38)
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

fname = paste("mcmcnbaout_",gammalen,"_",ndiv,".RData",sep='')
save.image(file=fname)



