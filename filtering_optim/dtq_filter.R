driftfun <- function(theta, y) 
{
  # computes drift as a matrix
  f = (theta[1])*(theta[2] - y)
  return(f)
}

difffun <- function(theta, y)
{
  # replicates diffusion as a long vector
  g = rep(theta[3], length(y))
  return(g)
}

difffun1 <- function(theta, y)
{
  # replicates diffusion term as a matrix
  dd1 = dim(y)[1]
  dd2 = dim(y)[2]
  vec = rep(theta[3], dd1)
  return(replicate(dd2, vec)) 
}

# source('integrandmat.R')
source('Dtheta.R')

dtq <- function(paramvec, h, k, bigm, littlet, data)
{
  # paramvec = {theta, x, sigeps2}
  numsteps = ceiling(littlet/h)
  xvec = c((-bigm):bigm)*k
  npts = length(xvec)
  
  A = integrandmat(xvec, xvec, h, driftfun, difffun, theta)
  D = Dtheta(xvec, xvec, h, driftfun, difffun1, theta)
  
  pdfmatrix = matrix(0,nrow=npts,ncol=(ncol(data)-1))
  
  ntheta = length(theta)
  qmattheta = list(NULL)
  for (i in c(1:ntheta)) qmattheta[[i]] = 0*pdfmatrix
  
  for (curcol in c(1:(ncol(data)-1)))
  {
    pdfmatrix[,curcol] = rowMeans(integrandmat(xvec, data[,curcol], h, driftfun, difffun, theta))
    initderivs = Dtheta(xvec, data[,curcol], h, driftfun, difffun1, theta)
    for (i in c(1:ntheta))
      qmattheta[[i]][,curcol] = rowMeans(initderivs[[i]])
  }
  startstep = 2
  
  for (curstep in c(startstep:(numsteps-1)))
  {
    pdfmatrix = k*A %*% pdfmatrix
    for (i in c(1:ntheta))
      qmattheta[[i]] = k*A %*% qmattheta[[i]] + k*D[[i]] %*% pdfmatrix
  }
  
  likelihood = matrix(0, nrow = nrow(data), ncol = (ncol(data)-1))
  gradient = list(NULL)
  for (i in c(1:ntheta)) gradient[[i]] = 0*likelihood
  
  for (curcol in c(2:ncol(data)))
  {
    # evaluate \Gamma vector across all samples at a particular time
    # this is a matrix because we're also evaluating at xvec
    # also evaluate the derivative of the \Gamma vector across all the samples and evaulated at the grid points
    
    gammamat = integrandmat(data[, curcol], xvec, h, driftfun, difffun, theta)
    likelihood[,(curcol-1)] = k*gammamat %*% pdfmatrix[,(curcol-1)]
    gdmat = Dtheta(data[,curcol], xvec, h, driftfun, difffun1, theta)
    for (i in c(1:ntheta))
      gradient[[i]][,(curcol-1)] = k*gdmat[[i]] %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta[[i]][,(curcol-1)]
  }
  
  return(list(lik = likelihood, grad = gradient))
}