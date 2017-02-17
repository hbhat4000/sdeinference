driftfun <- function(c0, y) 
{
  # computes drift as a matrix
  f = (c0[1])*(c0[2] - y)
  return(f)
}

difffun <- function(c0, y)
{
  # replicates diffusion as a long vector
  g = rep(c0[3], length(y))
  return(g)
}

difffun1 <- function(c0, y)
{
  # replicates diffusion term as a matrix
  dd1 = dim(y)[1]
  dd2 = dim(y)[2]
  vec = rep(c0[3], dd1)
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
    print(A[,2])
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
    print(qmattheta)
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


