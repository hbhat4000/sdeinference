driftfun <- function(c0,y) 
{
  return(0.5*(c0[1] - y))
}

difffun <- function(c0,y)
{
  return(rep(c0[2],length(y)))
}

source('integrandmat.R')

cdt <- function(c0,h,littlet,data)
{
    numsteps = ceiling(littlet/h)
    s = 0.75
    k = h^s
    yM = pi/k
    bigm = ceiling(yM/k)
    xvec = c((-bigm):bigm)*k
    npts = length(xvec)

    A = integrandmat(xvec,xvec,h,driftfun,difffun,c0)

    # pdf for first Dirac Delta initial condition, after one timestep
    newxvec = rep(xvec, dim(data)[2])
    thisdrift = driftfun(c0,data)
    meanvec1 = data+thisdrift*h
    meanvec = rep(meanvec1, each = npts)
    thisdiff = difffun(c0,data)
    sdvec1 = thisdiff*sqrt(h)
    sdvec = rep(sdvec1,each = npts)
    approxpdfsvec = dnorm(newxvec, mean=meanvec, sd=sdvec)
    pdfmatrix  = matrix(approxpdfsvec,nrow=npts,ncol=dim(data)[2])  


    # evolve all pdf's forward by correct # of time steps
    for (curstep in c(2:numsteps))
        pdfmatrix = k*A %*% pdfmatrix

    pdfmatrix[,(2:(ncol(pdfmatrix)))] = k*A %*% pdfmatrix[,(2:(ncol(pdfmatrix)))] 

    return(list(x=xvec,y=pdfmatrix))
}


