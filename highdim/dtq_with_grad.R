driftfun <- function(c0, y) 
{
  # computes drift as a matrix
  nc0 = length(c0) - 1
  if (class(y)=="numeric") y = as.matrix(y) 
  f = matrix(data = c0[1],nrow=nrow(y),ncol=ncol(y))
  for (i in c(2:nc0))
    f = f + c0[i]*y^(i-1)

  # f = c0[1]*y + c0[2]*y^3

  return(f)
}

difffun <- function(c0, y)
{
  # replicates diffusion as a long vector
  nc0 = length(c0)
  g = rep(c0[nc0], length(y))
  return(g)
}

difffun1 <- function(c0, y)
{
  # replicates diffusion term as a matrix
  dd1 = dim(y)[1]
  dd2 = dim(y)[2]
  nc0 = length(c0)
  vec = rep(c0[nc0], dd1)
  return(replicate(dd2, vec)) 
}

source('integrandmat.R')
source('Dtheta.R')

cdt <- function(c0, h, k, bigm, littlet, data)
{
    numsteps = ceiling(littlet/h)
    xvec = c((-bigm):bigm)*k
    npts = length(xvec)

    A = integrandmat(xvec, xvec, h, driftfun, difffun, c0)
    D = Dtheta(xvec, xvec, h, driftfun, difffun1, c0)
    
    pdfmatrix = matrix(0,nrow=npts,ncol=(ncol(data)-1))

    nc0 = length(c0)
    qmattheta = list(NULL)
    for (i in c(1:nc0)) qmattheta[[i]] = 0*pdfmatrix

    for (curcol in c(1:(ncol(data)-1)))
    {
        pdfmatrix[,curcol] = rowMeans(integrandmat(xvec, data[2:nrow(data),curcol], h, driftfun, difffun, c0))
        initderivs = Dtheta(xvec, data[2:nrow(data),curcol], h, driftfun, difffun1, c0)
        for (i in c(1:nc0))
            qmattheta[[i]][,curcol] = rowMeans(initderivs[[i]])
    }
    startstep = 2

    for (curstep in c(startstep:(numsteps-1)))
    {
        pdfmatrix = k*A %*% pdfmatrix
        for (i in c(1:nc0))
            qmattheta[[i]] = k*A %*% qmattheta[[i]] + k*D[[i]] %*% pdfmatrix
    }

    likelihood = matrix(0, nrow = (nrow(data) - 1), ncol = (ncol(data)-1))
    gradient = list(NULL)
    for (i in c(1:nc0)) gradient[[i]] = 0*likelihood

    for (curcol in c(2:ncol(data)))
    {
        # evaluate \Gamma vector across all samples at a particular time
        # this is a matrix because we're also evaluating at xvec
              
        gammamat = integrandmat(data[2:nrow(data), curcol], xvec, h, driftfun, difffun, c0)
        likelihood[,(curcol-1)] = k*gammamat %*% pdfmatrix[,(curcol-1)]
        gdmat = Dtheta(data[2:nrow(data),curcol], xvec, h, driftfun, difffun1, c0)
        for (i in c(1:nc0))
            gradient[[i]][,(curcol-1)] = k*gdmat[[i]] %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta[[i]][,(curcol-1)]
    }

    return(list(lik = likelihood, grad = gradient))
}


