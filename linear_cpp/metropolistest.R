rm(list = ls(all = TRUE))

# load the dtq package for finding the computing likelihood function
library('Rdtql')

# load data
load('fakedata.RData')
fd = xtraj

# define log prior
logprior <- function(z)
{
  return(dnorm(x = z, mean = 0, sd = 100, log = TRUE))
}

objfun <- function(theta, prior)
{
  numparam = length(theta)
  # probmat = array(0, dim = c(nrow(fd), ncol(fd) - 1, numparam + 1))
  probmat = Rdtql(method = "MCMC", thetavec = theta, h = myh, k = myk, M = mybigm, littlet = 1, init_data = fd)
  mylik = probmat[,,1]
  mylik[mylik < 0] = 0
  
  objective = sum(log(mylik)) + prior
  return(objective)
}

# create grid, compute densities on that grid
myh = 0.05
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.5))

# initial condition fakedata = c(1,4,0.5)
theta = c(1, 2, 1)
numparam = length(theta)

hh = 0.01
totsteps = 1000
thetamat = matrix(nrow=totsteps, ncol=numparam)
artrack = numeric(length=totsteps)

for (i in c(1:totsteps))
{
  # generate proposal 
  z = rnorm(n = numparam, mean = 0, sd = 0.25)
  theta_prop = theta + z
  
  oldlogprior = sum(logprior(prop))
  oldgradpost = objgradfun(theta, phi, mu, mass_sd, oldlogprior)
  
  phi_half = phi + (oldgradpost$gradient) * (hh / 2)
  theta_star = theta + phi_half * (hh / mass_sd)
  
  halflogprior = sum(logprior(phi_half, mu, mass_sd))
  propgradpost = objgradfun(theta_star, phi_half, mu, mass_sd, halflogprior)
  
  if (propgradpost$objective == -Inf)
    rho = 0
  else
  {
    phi_star = phi_half + (propgradpost$gradient) * (hh / 2)
    proplogprior = sum(logprior(phi_star, mu, mass_sd))
    rho = exp(propgradpost$objective + proplogprior - (oldgradpost$objective + oldlogprior))
  }
  
  # accept/reject step
  u = runif(n = 1)
  if (rho > u)
  {
    theta = theta_star
    oldgradpost = propgradpost
    print(paste("Accepted",paste("theta[",c(1:3),"]=",format(theta_star,digits=3,scientific=TRUE),collapse=', ',sep='')))
    artrack[i] = 1
  }
  else
  {
    print(paste("Rejected",paste("theta[",c(1:3),"]=",format(theta_star,digits=3,scientific=TRUE),collapse=', ',sep='')))
    artrack[i] = 0
  }
  thetamat[i,] = theta
}
myout = list(theta=thetamat,ar=artrack)
fname = paste('posteriorsamples_',myh,'.RData',sep='')
save(myout,file=fname)

