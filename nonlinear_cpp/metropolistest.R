rm(list = ls(all = TRUE))

# load all required parameters (initial theta, grid size, ...)
source('parameters.R')

# load the dtq package for finding the computing likelihood function
library('Rgdtq')

# load library for RMSE calculation
library('Metrics')

# load data as fd vector
fname = paste('fakedata_', datanum, '.RData', sep = '')
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
  mylik = Rdtq(thetavec = theta, h = gridh, k = gridk, M = gridM, littlet = 1, init_data = fd)
  mylik[mylik < 0] = 0
  
  objective = sum(log(mylik)) + prior
  return(objective)
}

# matrices which are used to store the computed parameters
thetamat = matrix(nrow = totsteps, ncol = numparam)
artrack = numeric(length = totsteps)
rmserror = numeric(length = totsteps)

ptm = proc.time()

# looping over
for (i in c(1:totsteps))
{
  oldlogprior = logprior(theta)
  oldlogpost = logpost(theta, oldlogprior)

  # generate proposal 
  z = rnorm(n = numparam, mean = prop_mu, sd = prop_sd)
  proptheta = theta + z
  
  proplogprior = logprior(proptheta)
  proplogpost = logpost(proptheta, proplogprior)
  
  rho = exp(proplogpost - oldlogpost)
  
  # accept/reject step
  u = runif(n = 1)
  if (rho > u)
  {
    theta = proptheta
    # oldlogprior = proplogprior
    # oldlogpost = proplogpost

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

myout = list(thetamat, artrack, rmserror, looptime)
fname = paste('MCMCsamples_', datanum, '.RData', sep = '')
save(myout, file = fname)


