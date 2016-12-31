f <- function(theta, x) 
{
  return((theta[1])*(theta[2] - x))
}

g <- function(theta, x)
{
  return(theta[3])
}

integrandmat <- function(x, y, h, f, g, theta)
{
  return(exp(-(x - y - f(theta,y)*h)^2/(2*g(theta,y)^2*h))/(abs(g(theta, y))*sqrt(2*pi*h)))
}

# front propagation from x_{i} to x_{i+1}
dtq_complete_front <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
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
  return(sum(log(approxpdf)))
}

# back propagation from x_{i} to x_{i+1}
dtq_complete_back <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
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
  return(sum(log(approxpdf)))
}

# one step Gaussian from x_{i} to z_{i1}: part1 = p(z_{i1} | x_{i})
# front propagation from z_{i1} to x_{i+1}: part2 = p(x_{i+1} | z_{i1})
dtq_firststep_front <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  part1 = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  
  part2 = k * (A %*% part1)
  
  for(i in c(2:numsteps-1))
    part2 = k * (A %*% part2)
  
  part2 = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta))) %*% part2
  
  # part1[part1 <= 2.2e-16] = 0
  part2[part2 <= 2.2e-16] = 0
  
  print(sum(log(part1)))
  print(sum(log(part2)))
  
  finalval = sum(log(part1)) * sum(log(part2))
  return(finalval)
}

# one step Gaussian from x_{i} to z_{i1}: part1 = p(z_{i1} | x_{i})
# back propagation from x_{i+1} to z_{i1}: part2 = p(x_{i+1} | z_{i1})
dtq_firststep_back <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  part1 = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  
  part2 = k * t(as.matrix(integrandmat(final, grid, h, f, g, theta)))
  
  # main loop
  for (i in c(2:numsteps-2))
    part2 = k * (part2 %*% A)
  
  part2 = k * (part2 %*% part1)
  
  # part1[part1 <= 2.2e-16] = 0
  part2[part2 <= 2.2e-16] = 0
  
  print(sum(log(part1)))
  print(sum(log(part2)))
  
  finalval = sum(log(part1)) * sum(log(part2))
  return(finalval)
}

# front propagation from x_{i} to z_{ij}: part1 = p(z_{ij} | x_{i})
# one step Gaussian from z_{ij} to z_{i,j+1}: part2 = p(z_{i,j+1} | z_{ij})
# front propagation from z_{i,j+1} to x_{i+1}: part3 = p(x_{i+1} | z_{i,j+1})
dtq_internal_front <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  

  # 
  # part1[part1 <= 2.2e-16] = 0
  # part2[part2 <= 2.2e-16] = 0
  # part3[part3 <= 2.2e-16] = 0
  # 
  # finalval = sum(log(part1)) * sum(log(part2) * sum(log(part3)))
  # return(finalval)
}

# front propagation from x_{i} to z_{ij}: part1 = p(z_{ij} | x_{i})
# one step Gaussian from z_{ij} to z_{i,j+1}: part2 = p(z_{i,j+1} | z_{ij})
# back propagation from x_{i+1} to z_{i,j+1}: part3 = p(x_{i+1} | z_{i,j+1})
dtq_internal_back <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
 
  # part1[part1 <= 2.2e-16] = 0
  # part2[part2 <= 2.2e-16] = 0
  # part3[part3 <= 2.2e-16] = 0
  # 
  # finalval = sum(log(part1)) * sum(log(part2) * sum(log(part3)))
  return(finalval)
}

# front propagation from x_{i} to z_{iF}: part1 = p(z_{iF} | x_{i})
# one step Gaussian from z_{iF} to x_{i+1}: part2 = p(x_{i+1} | z_{iF})
dtq_laststep_front <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)

  return(finalval)
}

# front propagation from x_{i} to z_{iF}: part1 = p(z_{iF} | x_{i})
# one step Gaussian back from x_{i+1} to z_{iF}: part2 = p(x_{i+1} | z_{iF})
dtq_laststep_back <- function(theta, h, k, M, deltat, init, final)
{
  numsteps = ceiling(deltat/h)
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  
  A = integrandmat(gridmat, t(gridmat), h, f, g, theta)
  
  part1 = k * (as.matrix(integrandmat(grid, init, h, f, g, theta)))
  
  for(i in c(2:numsteps-2))
    part1 = k * (A %*% part1)
  
  part2 = k * (as.matrix(integrandmat(final, grid, h, f, g, theta)))
    
  part1[part1 <= 2.2e-16] = 0
  part2[part2 <= 2.2e-16] = 0
  
  finalval = sum(log(part1)) * sum(log(part2))
  return(finalval)
}
