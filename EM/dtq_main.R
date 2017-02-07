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
  
  # \tau_{m}
  approxpdf = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))

  # A^{n-2} * \tau_{m}
  for (i in c(2:(numsteps-1)))
    approxpdf = k * (A %*% approxpdf)
  
  # \tau^T_{m+1} * A^{n-2} * \tau_{m}
  approxpdf = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta))) %*% approxpdf

  approxpdf[approxpdf <= 2.2e-16] = 0
  print(sum(log(approxpdf)))
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

# # one step Gaussian from x_{i} to z_{i1}: lambda = p(z_{i1} | x_{i})
# # front propagation from z_{i1} to x_{i+1}: gamma = p(x_{i+1} | z_{i1})
# dtq_firststep_front <- function(theta, h, k, M, numsteps, init, final) {
#   grid = c((-M):M)*k
#   gridmat = replicate(length(grid), grid)
  
#   A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
#   lambda = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  
#   # gamma = k * (A %*% lambda)
#   gamma = lambda
  
#   for(i in c(3:numsteps-1))
#     gamma = k * (A %*% gamma)
  
#   gamma = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta))) %*% gamma
  
#   # lambda[lambda <= 2.2e-16] = 0
#   # gamma[gamma <= 2.2e-16] = 0
  
#   print(c(sum(log(lambda)), sum(log(gamma))))
#   finalval = sum(log(lambda)) * sum(log(gamma))
#   return(finalval)
# }

# one step Gaussian from x_{i} to z_{i1}: lambda = p(z_{i1} | x_{i})
# back propagation from x_{i+1} to z_{i1}: gamma = p(x_{i+1} | z_{i1})
dtq_firststep_back <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  # lambda_i  
  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  
  # (1/k) \gamma_i A^{F-1}
  gamma = (1/k) * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  
  # main loop
  for (i in c(1:numsteps-1))
    gamma = k * (gamma %*% A)
  
  finalval = lambda % gamma
  logfinalval = sum(log(finalval))
  print(logfinalval)
  
  return(logfinalval)
}

# # front propagation from x_{i} to z_{ij}: lambda = p(z_{ij} | x_{i})
# # one step Gaussian from z_{ij} to z_{i,j+1}: gamma = p(z_{i,j+1} | z_{ij})
# # front propagation from z_{i,j+1} to x_{i+1}: part3 = p(x_{i+1} | z_{i,j+1})
# dtq_internal_front <- function(theta, h, k, M, numsteps, init, final) {
#   grid = c((-M):M)*k
#   gridmat = replicate(length(grid), grid)
  
#   A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
#   # return(finalval)
# }

# front propagation from x_{i} to z_{ij}: lambda = p(z_{ij} | x_{i})
# one step Gaussian from z_{ij} to z_{i,j+1}: gamma = p(z_{i,j+1} | z_{ij})
# back propagation from x_{i+1} to z_{i,j+1}: part3 = p(x_{i+1} | z_{i,j+1})
dtq_internal_back <- function(theta, h, k, M, numsteps, init, final, j) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  modifiedA = (1/k) * A

  lambda = as.matrix(integrandmat(grid, init, h, f, g, theta))
  gamma = (1/k) * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))

  for (i in c(1:j-1)) {
  	lambda = k * (A %*% lambda)
  }

  return(finalval)
}

# front propagation from x_{i} to z_{iF}: lambda = p(z_{iF} | x_{i})
# one step Gaussian from z_{iF} to x_{i+1}: gamma = p(x_{i+1} | z_{iF})
dtq_laststep_front <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)

  lambda = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  for(i in c(2:numsteps-2))
    lambda = k * (A %*% lambda)
  
  gamma = k * t(as.matrix(integrandmat(grid, init, h, f, g, theta))) %*% lambda
    
  print(c(sum(log(lambda)), sum(log(gamma))))
  finalval = sum(log(lambda)) * sum(log(gamma))
  return(finalval)
}

# front propagation from x_{i} to z_{iF}: lambda = p(z_{iF} | x_{i})
# one step Gaussian back from x_{i+1} to z_{iF}: gamma = p(x_{i+1} | z_{iF})
dtq_laststep_back <- function(theta, h, k, M, numsteps, init, final) {
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  lambda = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  for(i in c(2:numsteps-2))
    lambda = k * (A %*% lambda)
  
  gamma = k * (as.matrix(integrandmat(final, grid, h, f, g, theta)))
    
  # lambda[lambda <= 2.2e-16] = 0
  # gamma[gamma <= 2.2e-16] = 0
  
  print(c(sum(log(lambda)), sum(log(gamma))))
  finalval = sum(log(lambda)) * sum(log(gamma))
  return(finalval)
}
