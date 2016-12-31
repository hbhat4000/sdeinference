driftfun <- function(c0, y) 
{
  # computes drift as a matrix
  f = (c0[1])*y*(c0[2] - y^2)
  return(f)
}

difffun <- function(c0, y)
{
  # replicates diffusion as a long vector
  diffcoeff = exp(c0[3])
  g = rep(diffcoeff, length(y))
  return(g)
}

difffun1 <- function(c0, y)
{
  # replicates diffusion term as a matrix
  diffcoeff = exp(c0[3])
  dd1 = dim(y)[1]
  dd2 = dim(y)[2]
  vec = rep(diffcoeff, dd1)
  return(replicate(dd2, vec)) 
}

source('integrandmat.R')

dtq <- function(c0, h, k, bigm, littlet, data)
{
    numsteps = ceiling(littlet/h)
    xvec = c((-bigm):bigm)*k
    npts = length(xvec)

    A = integrandmat(xvec, xvec, h, driftfun, difffun, c0)
    pdfmatrix = matrix(0,nrow=npts,ncol=(ncol(data)-1))
    nc0 = length(c0)

    # initialization
    for (curcol in c(1:(ncol(data)-1)))
    {
        pdfmatrix[,curcol] = rowMeans(integrandmat(xvec, data[2:nrow(data),curcol], h, driftfun, difffun, c0))
    }
    startstep = 2

    # main loop
    for (curstep in c(startstep:(numsteps-1)))
        pdfmatrix = k*A %*% pdfmatrix

    likelihood = matrix(0, nrow = (nrow(data) - 1), ncol = (ncol(data)-1))

    # the next step is instead of interpolation
    for (curcol in c(2:ncol(data)))
    {
        gammamat = integrandmat(data[2:nrow(data), curcol], xvec, h, driftfun, difffun, c0)
        likelihood[,(curcol-1)] = k*gammamat %*% pdfmatrix[,(curcol-1)]
    }
    return(likelihood)
}


