drift <- function(theta, y) 
{
  # computes drift as a matrix
  f = (theta[1])*(theta[2] - y)
  return(f)
}

diffusion <- function(theta, y)
{
  # fills y by theta[3]
  g = theta[3] + y*0
  return(g)
}

source('Dtheta.R')
source('integrandmat.R')

dtq <- function(theta, h, k, bigm, littlet, data)
{
  numsteps = ceiling(littlet/h)
  xvec = c((-bigm):bigm)*k
  npts = length(xvec)
  
  A = integrandmat(xvec, xvec, h, theta)
  D = Dtheta(xvec, xvec, h, theta)
  
  pdfmatrix = matrix(0, nrow = npts, ncol = (ncol(data)-1))
  
  ntheta = length(theta)
  qmattheta = list(NULL)
  for (i in c(1:ntheta)) qmattheta[[i]] = 0*pdfmatrix
  
  for (curcol in c(1:(ncol(data)-1)))
  {
    pdfmatrix[,curcol] = rowMeans(integrandmat(xvec, data[,curcol], h, theta))
    initderivs = Dtheta(xvec, data[,curcol], h, theta)
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
    
    gammamat = integrandmat(data[, curcol], xvec, h, theta)
    likelihood[,(curcol-1)] = k*gammamat %*% pdfmatrix[,(curcol-1)]
    gdmat = Dtheta(data[,curcol], xvec, h, theta)
    for (i in c(1:ntheta))
      gradient[[i]][,(curcol-1)] = k*gdmat[[i]] %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta[[i]][,(curcol-1)]
  }
  
  return(list(lik = likelihood, grad = gradient))
}