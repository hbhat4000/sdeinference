f <- function(theta, x, returngrad=FALSE) {
  fval = (theta[1])*(theta[2] - x)
  if (returngrad==FALSE)
    return(fval)
  else
  {
    if (class(x)=="numeric")
    {
      dfdtheta = matrix(0,nrow=length(x),ncol=length(theta))
      dfdtheta[,1] = fval/theta[1]
      dfdtheta[,2] = theta[1]
      dfdtheta[,3] = 0
    }
    else if (class(x)=="matrix")
    {
      dfdtheta = array(0,c(dim(x),length(theta)))
      dfdtheta[,,1] = fval/theta[1]
      dfdtheta[,,2] = theta[1]
      dfdtheta[,,3] = 0
    }
    return(list(val=fval,grad=dfdtheta))
  } 
}

g <- function(theta, x, returngrad=FALSE) {
  gval = theta[3] # warning this will not have the same size as x
  if (returngrad==FALSE)
    return(gval)
  else
  {
    dgdtheta = c(0,0,1)
    #dgdtheta[,1] = 0
    #dgdtheta[,2] = 0
    #dgdtheta[,3] = 1
    return(list(val=gval,grad=dgdtheta))
  }
}

integrandmat <- function(x, y, h, f, g, theta) {
  return(exp(-(x - y - f(theta, y)*h)^2/(2*g(theta, y)^2*h))/(abs(g(theta, y))*sqrt(2*pi*h)))
}

# front propagation from x_{i} to x_{i+1}
dtq_complete_front <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  savepdf = matrix(0,ncol=(numsteps-1),nrow=(2*M+1))
  
  # \tau_{m}
  if(numsteps >= 1) {
    # this is actually lambda, which has no k
    approxpdf = (as.matrix(integrandmat(grid, init, h, f, g, theta)))
    savepdf[,1] = approxpdf
  }
  
  # A^{n-2} * \tau_{m}
  # note that k * A = K in the notes
  if(numsteps >= 3) {
    for (i in c(2:(numsteps-1))) {
      approxpdf = k * (A %*% approxpdf)
      savepdf[,i] = approxpdf
    }
  }
  
  # \tau^T_{m+1} * A^{n-2} * \tau_{m}
  gamma = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  if (numsteps >= 2) {
    approxpdf = gamma %*% approxpdf
  }
  
  # approxpdf[approxpdf <= 2.2e-16] = 0
  #print(sum(log(approxpdf)))
  #savepdf[savepdf <= 2.2e-16] = 0
  
  #print(savepdf)
  #plot(savepdf)
  
  return(list(lik=as.numeric(approxpdf),pdf=savepdf))
}

# back propagation from x_{i} to x_{i+1}
dtq_complete_back <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  # \tau^T_{m+1}
  # this is actually gamma
  approxpdf = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  
  # \tau^T_{m+1} * A^{n-2}
  # again, k * A corresponds to K in the notes
  for (i in c(2:(numsteps-1)))
    approxpdf = k * (approxpdf %*% A)
  
  # \tau^T_{m+1} * A^{n-2} * \tau_{m}
  # final multiplication by lambda, note that lambda has no k
  approxpdf = approxpdf %*% (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  
  # approxpdf[approxpdf <= 2.2e-16] = 0
  return(as.numeric(approxpdf))
}

# one step Gaussian from x_{i} to z_{i1}: lambda = p(z_{i1} | x_{i})
# back propagation from x_{i+1} to z_{i1}: gamma = p(x_{i+1} | z_{i1})
# note: this function returns a pdf
dtq_firststep <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  # \lambda_i  
  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  
  # \gamma_i
  gamma = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  
  # (1/k) \gamma_i A^{F-1}
  # numsteps = F+1, so F-1 = numsteps-2
  for (i in c(1:(numsteps-2)))
    gamma = k * (gamma %*% A)
  
  finalpdf = (1/k)*hadamard.prod(t(gamma), lambda)
  logfinalpdf = log(finalpdf)
  # print(logfinalval)
  
  return(finalpdf)
}

# front propagation from x_{i} to z_{ij}: lambda = p(z_{ij} | x_{i})
# one step Gaussian from z_{ij} to z_{i,j+1}: gamma = p(z_{i,j+1} | z_{ij})
# back propagation from x_{i+1} to z_{i,j+1}: part3 = p(x_{i+1} | z_{i,j+1})
# note: this function returns a two-dimensional pdf as a matrix
dtq_internal <- function(theta, h, k, M, numsteps, init, final, j) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  # modifiedA = (1/k) * A

  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  gamma = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))

  if (j >= 2)
  {
    for (i in c(1:(j-1))) {
      lambda = k * (A %*% lambda)
    }
  }
    
  # numsteps = F+1, so F = numsteps-1
  if (j <= (numsteps-3))
  {
    for (i in c(1:(numsteps-j-2))) {
      gamma = k * (gamma %*% A)
    }
  }
    
  kronpdt = t(kronecker(gamma, lambda))
  finalpdf2 = (1/k)*hadamard.prod(kronpdt, A)
  return(finalpdf2)
}

# front propagation from x_{i} to z_{iF}: lambda = p(z_{iF} | x_{i})
# one step Gaussian back from x_{i+1} to z_{iF}: gamma = p(x_{i+1} | z_{iF})
# note: this function returns a pdf
dtq_laststep <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  for(i in c(1:(numsteps-2)))
    lambda = k * (A %*% lambda)
  
  gamma = k * (as.matrix(integrandmat(final, grid, h, f, g, theta)))

  finalpdf = (1/k)*hadamard.prod(gamma, lambda)
  logfinalpdf = log(finalpdf)
  # print(logfinalval)
  return(finalpdf)
}

# we want this function to return three things:
# 1) p(z_{i1} | x, \theta^k)
# 2) p(z_{i,j+1}, z_{ij} | x, \theta^k)
# 3) p(z_{iF} | x, \theta^k)
dtq_all <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)

  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  gamma = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  
  # the list index is 1 + the power of K that we have multiplied by
  lambdalist = list(NULL)
  lambdalist[[1]] = lambda
  gammalist = list(NULL)
  gammalist[[1]] = gamma
  
  # numsteps-2 = F-1, which is the highest power of K that we have
  for (j in c(1:(numsteps-2)))
  {
    lambdalist[[j+1]] = (k*A) %*% lambdalist[[j]]
    gammalist[[j+1]] = gammalist[[j]] %*% (k*A)
  }
  
  # compute the complete likelihood, p(x_{i+1} | x_i, \theta)
  complete = as.numeric(gammalist[[1]] %*% lambdalist[[numsteps-1]])
  
  # compute first and last 1d pdf's
  first = (1/k)*hadamard.prod(t(gammalist[[numsteps-1]]),lambdalist[[1]])/complete
  last = (1/k)*hadamard.prod(t(gammalist[[1]]),lambdalist[[numsteps-1]])/complete
  
  # compute all intermediate 2d pdf's
  pdf2dlist = list(NULL)
  for (j in c(1:(numsteps-2))) 
  {
    pdf2dlist[[j]] = hadamard.prod(t(kronecker(gammalist[[(numsteps-1-j)]],lambdalist[[j]])),A)/(k*complete)
  }    
  return(list(complete=complete,first=first,last=last,pdf2d=pdf2dlist))
}
