rm(list = ls(all = TRUE))

# load data
load('fakedata.RData')
fd = xtraj

# load necessary functions
source('dtq_filter.R')

objgradfun <- function(theta0)
{
  probmat = dtq(theta0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)
  mylik = probmat$lik
  mylik[mylik < 0] = 0
  objective = -sum(log(mylik))
  ntheta0 = length(theta0)
  gradient = numeric(length=ntheta0)
  for (i in c(1:ntheta0))
    gradient[i] = -sum(probmat$grad[[i]] / probmat$lik)
  
  return(list("objective" = objective, "gradient" = gradient))
}

# create grid, compute densities on that grid
myh = 0.05
myk = myh^(0.75)
mybigm = ceiling(pi/(myk^1.5))

theta = c(2, 2, 1)

for(i in c(1:length(xtraj)))
  # optimize using theta, x_n, sigma_epsilon2
  paramvec = c(theta, xtraj[1], sigeps2)

library('nloptr')

# print(objgradfun(theta))

res <- nloptr(x0 = paramvec, eval_f = objgradfun, lb = c(0.1, 0, 0.1), ub = c(10, 10, 2), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4))