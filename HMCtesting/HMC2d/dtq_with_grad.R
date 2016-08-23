source('driftdiff.R')
source('integrandmat.R')
# inplace of integrandmat
library('Rdtq2d')
source('Dtheta.R')

dtq <- function(thetavec, h, k, bigm, littlet, data)
{
  x1 = data[[2]]
  x2 = data[[3]]

  numsteps = ceiling(littlet/h)
  xvec = c((-bigm):bigm)*k
  npts = length(xvec)

  #A = integrandmat(xvec, xvec, h, driftfun, difffun, thetavec)
  A = Rdtq2d(thetavec, x1, x2, h, numsteps, k, bigm)
  D = Dtheta(xvec, xvec, h, driftfun, difffun1, thetavec)
  
  pdfmatrix = matrix(0, nrow = npts, ncol = (ncol(data)-1))

  thetavec = length(thetavec)
  qmattheta = list(NULL)

  for (i in c(1:thetavec)) 
    qmattheta[[i]] = 0*pdfmatrix

  for (curcol in c(1:(ncol(data)-1)))
  {
    pdfmatrix[,curcol] = rowMeans(integrandmat(xvec, data[2:nrow(data),curcol], h, driftfun, difffun, thetavec))
    initderivs = Dtheta(xvec, data[2:nrow(data),curcol], h, driftfun, difffun1, thetavec)
    for (i in c(1:thetavec))
      qmattheta[[i]][,curcol] = rowMeans(initderivs[[i]])
  }
  startstep = 2

  for (curstep in c(startstep:(numsteps-1)))
  {
    pdfmatrix = k*A %*% pdfmatrix
    for (i in c(1:thetavec))
      qmattheta[[i]] = k*A %*% qmattheta[[i]] + k*D[[i]] %*% pdfmatrix
  }

  likelihood = matrix(0, nrow = (nrow(data) - 1), ncol = (ncol(data) - 1))
  gradient = list(NULL)

  for (i in c(1:thetavec)) 
    gradient[[i]] = 0*likelihood

  for (curcol in c(2:ncol(data)))
  {
    # evaluate \Gamma vector across all samples at a particular time
    # this is a matrix because we're also evaluating at xvec
          
    gammamat = integrandmat(data[2:nrow(data), curcol], xvec, h, driftfun, difffun, thetavec)
    likelihood[,(curcol-1)] = k*gammamat %*% pdfmatrix[,(curcol-1)]
    gdmat = Dtheta(data[2:nrow(data),curcol], xvec, h, driftfun, difffun1, thetavec)
    for (i in c(1:thetavec))
      gradient[[i]][,(curcol-1)] = k*gdmat[[i]] %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta[[i]][,(curcol-1)]
  }

  return(list(lik = likelihood, grad = gradient))
}


