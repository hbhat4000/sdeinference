rm(list = ls(all = TRUE))

# load data
load('fakedata.RData')
fd = ytraj

# load necessary functions
source('dtq_filter.R')

# define log prior which computes the sum of log normals (prior)
logprior <- function(paramvec)
{
  # paramvec = (theta[1], theta[2], theta[3], x, sigeps2)
  logpriortheta = dnorm(x = paramvec[1], mean = 0.5, sd = 1, log = TRUE) + dnorm(x = paramvec[2], mean = 2, sd = 10, log = TRUE)
  
  # y_0 = x_0 + \epsilon_0 => x_0 ~ N(y_0, \sigma^2_{\epsilon})
  logpriorx0 = dnorm(x = paramvec[4], mean = fd[1], sd = paramvec[5], log = TRUE)
  logpriorsigeps2 = dexp(x = paramvec[5], rate = 1, log = TRUE)

  return(logpriortheta + logpriorx0 + logpriorsigeps2)
}

logGaussian <- function(paramvec, xvec, yvec)
{
  # TODO: figure out which y and x to evaluate the gaussian at
  likGauss = dnorm(x = yvec, mean = xvec, sd = paramvec[5], log = TRUE)
  return(sum(likGauss))
}

Gaussianderivsigeps2 <- function(xvec, yvec, sigeps2)
{
  val = -log(sqrt(2*pi*sigeps2))/sigeps2 - (yvec - xvec)^2/(2*(sigeps2)^2)
  return(sum(val))
}

gammavecderiv <- function(xval, grid)
{
  part1 = (xval - grid - drift(grid, theta)*h)
  part2 = abs(diffusion(grid, theta))^2*h
  # G = exp(-(part1)^2/(2*part2))/(sqrt(2*pi*part2))
  Gderiv = -(part1)/(part2) #*G
  return(Gderiv)
}

# logposterior = log likelihood (computed by Gaussians) + log likelihood (computed by DTQ) + log prior 
# logposterior = log(y_j | x_j, \theta) + sum(log(DTQ output/mylik)) + logprior(x_0) + logprior(theta) + logprior(sigeps2)
# likelihood is computed using the DTQ method
logposterior <- function(paramvec)
{
  theta = paramvec[1:3]
  x = paramvec[4]
  sigeps2 = paramvec[5]
  grid = c((-mybigm):mybigm)*myk

  probmat = dtq(theta, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)

  DTQlik = probmat$lik
  DTQlik[DTQlik < 0] = 0

  priorval = logprior(paramvec)
  likval = logGaussian(paramvec)
  objective = sum(log(DTQlik)) + priorval + likval

  ntheta = length(theta)
  nparam = length(paramvec)
  gradient = numeric(length = nparam)

# TODO: figure out which y to use
  deriv1 = -(theta[1] - 0.5)
  deriv2 = -(theta[2] - 2)/100
  deriv3 = 0
  deriv4 = (y - x)/(sigeps2) + gammavecderiv(x, grid)
  deriv5 = Gaussianderivsigeps2(xvec, yvec, sigeps2) - 1

  # add the theta, x and sigeps2 derivs
  for (i in c(1:ntheta))
    gradient[i] = sum(probmat$grad[[i]] / probmat$lik)
  
  gradient[1] = gradient[1] + deriv1
  gradient[2] = gradient[2] + deriv2
  gradient[3] = gradient[3] + deriv3
  gradient[4] = deriv4
  gradient[5] = deriv5
  return(list("objective" = -objective, "gradient" = -gradient))
}

# create grid, compute densities on that grid
myh = 0.05
myk = myh^(0.75)
mybigm = ceiling(pi/(myk^1.5))

# optimize using theta, x_n, sigma_epsilon2
# paramvec = c(theta, ytraj[1], sigeps2)
# BigMine paper, initial theta = c(1.0,0.1,0.25), sigeps = 1.0
paramvec = (2, 2, 1, fd[1], 1)

library('nloptr')
res <- nloptr(x0 = paramvec, eval_f = logposterior, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4))