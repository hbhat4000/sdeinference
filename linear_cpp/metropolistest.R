rm(list = ls(all = TRUE))

# load all required parameters (initial theta, grid size, ...)
source('parameters.R')

# load the dtq package for finding the computing likelihood function
library('Rgdtq')

# load library for RMSE calculation
library('Metrics')

# load data as fd vector
fname = paste('fakedata_', fakedatanum, '.RData', sep = '')
load(fname)
fd = xtraj

# define log prior which is summed up before returning
logprior <- function(z)
{
  return(sum(dnorm(x = z, mean = prior_mu, sd = prior_sd, log = TRUE)))
}

# define log posterior as summation of log likelihood and log prior (computed in the function logprior).
# log likelihood is computed by calling the C++ function Rdtq which 
# returns the likelihood function computed by DTQ method
logpost <- function(theta, prior)
{
  mylik = Rdtql(thetavec = theta, h = gridh, k = gridk, M = gridM, littlet = 1, init_data = fd)
  mylik[mylik < 0] = 0
  
  objective = sum(log(mylik)) + prior
  return(objective)
}

# matrices which are used to store the computed parameters
thetamat = matrix(nrow = totsteps, ncol = numparam)
thetagenerated = matrix(nrow = totsteps, ncol = numparam)
artrack = numeric(length = totsteps)
rmserror = numeric(length = totsteps)

ptm = proc.time()

# looping over
oldlogprior = logprior(theta)
oldlogpost = logpost(theta, oldlogprior)

for (i in c(1:totsteps))
{
  # generate proposal 
  z = rnorm(n = numparam, mean = prop_mu, sd = prop_sd)
  proptheta = theta + z
  
  proplogprior = logprior(proptheta)
  proplogpost = logpost(proptheta, proplogprior)
  
  rho = exp(proplogpost - oldlogpost)
  
  thetagenerated[i,] = proptheta
  # accept/reject step
  u = runif(n = 1)
  if (rho > u)
  {
    theta = proptheta
    oldlogprior = proplogprior
    oldlogpost = proplogpost

    print(paste("Accepted step", i, ": ", paste("theta[", c(1:numparam), "]=", format(proptheta, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
    artrack[i] = 1
  }
  else
  {
    print(paste("Rejected step", i, ": ", paste("theta[", c(1:numparam), "]=", format(proptheta, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
    artrack[i] = 0
  }

  # Assigning values irrespective of accept/reject step
  thetamat[i,] = theta
  rmserror[i] = rmse(actualtheta, theta)
}

looptime = proc.time() - ptm

acceptance = 1 - mean(duplicated(thetamat[-(1:burnin),]))

par(mfrow = c(2,3))

hist(thetamat[-(1:burnin),1], nclass = 30, main = "Posterior of theta1", xlab = "True value = red line")
abline(v = mean(thetamat[-(1:burnin),1]))
abline(v = actualtheta[1], col = "red")

hist(thetamat[-(1:burnin),2], nclass = 30, main = "Posterior of theta2", xlab = "True value = red line")
abline(v = mean(thetamat[-(1:burnin),2]))
abline(v = actualtheta[2], col = "red")

hist(thetamat[-(1:burnin),3], nclass = 30, main = "Posterior of theta3", xlab = "True value = red line")
abline(v = mean(thetamat[-(1:burnin),3]))
abline(v = actualtheta[3], col = "red")

plot(thetamat[-(1:burnin),1], type = "l", xlab = "True value = red line", main = "MCMC values of theta1")
abline(h = actualtheta[1], col = "red")

plot(thetamat[-(1:burnin),2], type = "l", xlab = "True value = red line", main = "MCMC values of theta2")
abline(h = actualtheta[2], col = "red")

plot(thetamat[-(1:burnin),3], type = "l", xlab = "True value = red line", main = "MCMC values of theta3")
abline(h = actualtheta[3], col = "red")

# for comparison:
# summary(lm(y~x))

myout = list(gridh, gridk, gridM, actualtheta, burnin, totsteps, prop_mu, prop_sd, prior_mu, prior_sd, thetamat, artrack, rmserror, looptime)
fname = paste('MCMCsamples_', mcmcdatanum, '.RData', sep = '')
save(myout, file = fname)


