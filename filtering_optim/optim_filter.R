rm(list = ls(all = TRUE))

# load data
load('fakedata.RData')
fd = xtraj

# load necessary functions
source('dtq_filter.R')

# define log prior which computes the sum of log normals (prior)
logprior <- function(paramvec)
{
  logpriortheta = dnorm(x = theta[1], mean = 0.5, sd = 1, log = TRUE) + dnorm(x = theta[2], mean = 2, sd = 10, log = TRUE)
  # y_0 = x_0 + \epsilon_0 => x_0 ~ N(y_0, \sigma^2_{\epsilon})
  logpriorsigeps2 = dexp(x = sigeps2, rate = 1, log = TRUE)
  logpriorx0 = dnorm(x = , mean = )

  return(logpriortheta + logpriorx0 + logpriorsigeps2)
}

# logposterior = log likelihood + log prior 
# logposterior = log(y_j | x_j, \theta) + sum(log(DTQ output/mylik)) + logprior(x_0) + logprior(theta) + logprior(sigeps2)
# likelihood is computed using the DTQ method
logposterior <- function(paramvec)
{
  probmat = dtq(paramvec, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)
  mylik = probmat$lik
  mylik[mylik < 0] = 0
  objective = sum(log(mylik))
  nparam = length(paramvec)
  gradient = numeric(length = nparam)
  for (i in c(1:nparam))
    gradient[i] = sum(probmat$grad[[i]] / probmat$lik)
  
  return(list("objective" = -objective, "gradient" = -gradient))
}

# create grid, compute densities on that grid
myh = 0.05
myk = myh^(0.75)
mybigm = ceiling(pi/(myk^1.5))

theta = c(2, 2, 1)

# optimize using theta, x_n, sigma_epsilon2
# paramvec = c(theta, xtraj[1], sigeps2)
paramvec = list("theta" = theta, "xval" = xtraj, "sigeps2" = sigeps2)
library('nloptr')
res <- nloptr(x0 = paramvec, eval_f = objgradfun, lb = c(0.1, 0, 0.1), ub = c(10, 10, 2), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4))