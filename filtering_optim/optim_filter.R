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
  logpriorx0 = dnorm(x = paramvec[4], mean = 0, sd = 1, log = TRUE)
  logpriorsigeps2 = dexp(x = paramvec[5], rate = 1, log = TRUE)

  return(logpriortheta + logpriorx0 + logpriorsigeps2)
}

logGaussian <- function(paramvec, xvec, yvec)
{
  # TODO: figure out which y and x to evaluate the gaussian at
  likGauss = dnorm(x = yvec, mean = xvec, sd = sqrt(paramvec[5]), log = TRUE)
  return(sum(likGauss))
}

Gaussianderivsigeps2 <- function(xvec, yvec, sigeps2)
{
  # val = -log(sqrt(2*pi*sigeps2))/sigeps2 + (yvec - xvec)^2/(2*(sigeps2)^2)
  val = -1/(2*sigeps2) + (yvec - xvec)^2/(2*(sigeps2)^2)
  return(sum(val))
}

gammavecderiv <- function(xval, grid, theta, h)
{
  part1 = (xval - grid - drift(theta,grid)*h)
  part2 = abs(diffusion(theta,grid))^2*h
  G = exp(-(part1)^2/(2*part2))/(sqrt(2*pi*part2))
  Gderiv = -(part1)/(part2)*G
  return(Gderiv)
}

# logposterior = log likelihood (computed by Gaussians) + log likelihood (computed by DTQ) + log prior 
# logposterior = log(y_j | x_j, \theta) + sum(log(DTQ output/mylik)) + logprior(x_0) + logprior(theta) + logprior(sigeps2)
# likelihood is computed using the DTQ method
logposterior <- function(paramvec,thisx,thisy)
{
  theta = paramvec[1:3]
  x = paramvec[4]
  sigeps2 = paramvec[5]
  grid = c((-mybigm):mybigm)*myk

  # create observation series
  thisy1 = matrix(thisy[-1,],nrow=1)

  # number of parameters, initialization of gradient
  ntheta = length(theta)
  nparam = length(paramvec)
  gradient = numeric(length = nparam)

  # latest y value, corresponding to same time as "x" above
  lasty = thisy1[1,ncol(thisy1)] 

  objective = 0
  if (ncol(thisy1) >= 2)
  {
    # augment state series by last state (the one we're solving for!)
    thisx1 = matrix(thisx[-1,],nrow=1)
    thisx1 = cbind(thisx1,x)

    # time step of the data
    deltat = thisy[1,2] - thisy[1,1]  # assuming equispaced
    probmat = dtq(theta, h = myh, k = myk, bigm = mybigm, littlet = deltat, data = thisx1)

    DTQlik = probmat$lik
    DTQlik[DTQlik < 0] = 0

    objective = objective + sum(log(DTQlik))

    # add the theta, x and sigeps2 derivs
    for (i in c(1:ntheta))
      gradient[i] = sum(probmat$grad[[i]] / probmat$lik)

    # derivative of log p(x_{L+1} | x_L, \theta) w.r.t. x_{L+1}
    gammapiece = myk * (gammavecderiv(x, grid, theta, myh) %*% probmat$lv)
    gammapiece = gammapiece/probmat$lik[ncol(thisx1)-1]
    gradient[4] = gammapiece
  }
  else
  {
    # state series consists only of last state
    thisx1 = matrix(x,nrow=1,ncol=1)

    # if we are indeed solving for x0, then we need 
    # the derivative of p(x0) w.r.t. x0, instead of gammapiece above
    gradient[4] = -x
  }
  priorval = logprior(paramvec)
  likval = logGaussian(paramvec, thisx1, thisy1)
  objective = objective + priorval + likval


# TODO: figure out which y to use
  deriv1 = -(theta[1] - 0.5)
  deriv2 = -(theta[2] - 2)/100
  deriv3 = 0
  deriv4 = (lasty - x)/(sigeps2) 
  deriv5 = Gaussianderivsigeps2(thisx1, thisy1, sigeps2) - 1

  gradient[1] = gradient[1] + deriv1
  gradient[2] = gradient[2] + deriv2
  gradient[3] = gradient[3] + deriv3
  gradient[4] = gradient[4] + deriv4
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

library('nloptr')

mylb = c(-Inf,-Inf,sqrt(myk),-Inf,1e-4)

# res <- nloptr(x0 = paramvec, eval_f = logposterior, lb = mylb, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4),thisx=fd[,-ncol(fd)],thisy=fd)

# myfd = matrix(fd[,1],nrow=2)
# res <- nloptr(x0 = paramvec, eval_f = logposterior, lb = mylb, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4),thisx=NULL,thisy=myfd)

# iteratively run filtering+inference procedure
myx = NULL
paramvec = c(2, 2, 1, 0, 1)
for (j in c(1:ncol(fd)))
{
  myfd = matrix(fd[,1:j],nrow=2)
  paramvec[4] = myfd[2,j]
  res <- nloptr(x0 = paramvec, eval_f = logposterior, lb = mylb, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4),thisx=myx,thisy=myfd)
  if (is.null(myx))
  {
    myx = matrix(0,nrow=2,ncol=1)
    myx[1,1] = myfd[1,1]
    myx[2,1] = res$solution[4]
  }
  else
  {
    myx = cbind(myx,c(myfd[1,j],res$solution[4]))
  }
  paramvec = res$solution
}






