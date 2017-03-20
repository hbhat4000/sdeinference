f <- function(theta, x) {
  return((theta[1])*(theta[2] - x))
}

g <- function(theta, x) {
  return(theta[3])
}

integrandmat <- function(x, y, h, f, g, theta) {
  return(exp(-(x - y - f(theta, y)*h)^2/(2*g(theta, y)^2*h))/(abs(g(theta, y))*sqrt(2*pi*h)))
}

# front propagation from x_{i} to x_{i+1}
dtq_complete_front <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  savepdf = numeric(numsteps)
  
  # \tau_{m}
  if(numsteps >= 1) {
    approxpdf = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))
    savepdf[1] = approxpdf[M + 1]
  }
  
  # A^{n-2} * \tau_{m}
  if(numsteps >= 3) {
    for (i in c(2:(numsteps-1))) {
      approxpdf = k * (A %*% approxpdf)
      savepdf[i] = approxpdf[M + 1]
    }
  }
  
  # \tau^T_{m+1} * A^{n-2} * \tau_{m}
  if (numsteps >= 2) {
    approxpdf = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta))) %*% approxpdf
    savepdf[numsteps] = approxpdf[M + 1]  
  }
  
  # print(approxpdf)
  approxpdf[approxpdf <= 2.2e-16] = 0
  print(sum(log(approxpdf)))
  
  plot(savepdf)
  
  return(sum(log(approxpdf)))
}

# back propagation from x_{i} to x_{i+1}
dtq_complete_back <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  # \tau^T_{m+1}
  approxpdf = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  
  # \tau^T_{m+1} * A^{n-2}
  for (i in c(2:(numsteps-1)))
    approxpdf = k * (approxpdf %*% A)
  
  # \tau^T_{m+1} * A^{n-2} * \tau_{m}
  approxpdf = k * approxpdf %*% (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  
  approxpdf[approxpdf <= 2.2e-16] = 0
  print(sum(log(approxpdf)))
  return(sum(log(approxpdf)))
}

# one step Gaussian from x_{i} to z_{i1}: lambda = p(z_{i1} | x_{i})
# back propagation from x_{i+1} to z_{i1}: gamma = p(x_{i+1} | z_{i1})
dtq_firststep <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  # \lambda_i  
  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  
  # \gamma_i
  gamma = (1/k) * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  
  # (1/k) \gamma_i A^{F-1}
  for (i in c(1:numsteps-1))
    gamma = k * (gamma %*% A)
  
  finalval = hadamard.prod(t(gamma), lambda)
  logfinalval = sum(log(finalval))
  # print(logfinalval)
  
  return(logfinalval)
}

# front propagation from x_{i} to z_{ij}: lambda = p(z_{ij} | x_{i})
# one step Gaussian from z_{ij} to z_{i,j+1}: gamma = p(z_{i,j+1} | z_{ij})
# back propagation from x_{i+1} to z_{i,j+1}: part3 = p(x_{i+1} | z_{i,j+1})
dtq_internal <- function(theta, h, k, M, numsteps, init, final, j) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  modifiedA = (1/k) * A

  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  gamma = (1/k) * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))

  for (i in c(1:j-1)) {
  	lambda = k * (A %*% lambda)
  }

  for (i in c(1:numsteps-(j+1))) {
  	gamma = k * (gamma %*% A)
  }

  kronpdt = kronecker(gamma, lambda)
  
  finalval = hadamard.prod(kronpdt, modifiedA)
  logfinalval = sum(log(finalval))
  # print(c(j,logfinalval))
  return(logfinalval)
}

# front propagation from x_{i} to z_{iF}: lambda = p(z_{iF} | x_{i})
# one step Gaussian back from x_{i+1} to z_{iF}: gamma = p(x_{i+1} | z_{iF})
dtq_laststep <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  for(i in c(1:numsteps-1))
    lambda = k * (A %*% lambda)
  
  gamma = (1/k) * (as.matrix(integrandmat(final, grid, h, f, g, theta)))

  finalval = hadamard.prod(gamma, lambda)
  logfinalval = sum(log(finalval))
  # print(logfinalval)
  return(logfinalval)
}
