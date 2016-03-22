driftfun <- function(c0,y) 
{
  return(0.5*(c0[1] - y))
}

difffun <- function(c0,y)
{
  return(rep(c0[2],length(y)))
}

difffun1 <- function(c0,y)
{
  dd1 = dim(y)[1]
  dd2 = dim(y)[2]
  vec = rep(c0[2],dd1)
  return(replicate(dd2, vec)) 
}


source('integrandmat.R')
# source('diffgandf.R')
source('Dtheta1.R')
source('Dtheta2.R')
source('diffdiffustheta1.R')
source('diffdrifttheta1.R')
source('diffdiffustheta2.R')
source('diffdrifttheta2.R')

cdt <- function(c0,h,littlet,pdfmatrix0,data)
{
    numsteps = ceiling(littlet/h)
    s = 0.75
    k = h^s
    yM = pi/(k^0.5)
    bigm = ceiling(yM/k)
    xvec = c((-bigm):bigm)*k
    npts = length(xvec)

    A = integrandmat(xvec,xvec,h,driftfun,difffun,c0)
    D1 = Dtheta1(xvec,xvec,h,driftfun,difffun1,c0) # Derivative of A matrix w.r.t. theta1
    D2 = Dtheta2(xvec,xvec,h,driftfun,difffun1,c0) # Derivative of A matrix w.r.t. theta2

    # deravatives of the pdf w.r.t theta1 and theta2
    qmattheta1 = 0*pdfmatrix0
    qmattheta2 = 0*pdfmatrix0

    # evolve all pdf's forward by correct # of time steps
    pdfmatrix = pdfmatrix0
    likelihood = matrix(0,nrow=(nrow(data)-1),ncol=(ncol(data)-1))
    d1mat = matrix(0,nrow=(nrow(data)-1),ncol=(ncol(data)-1))
    d2mat = matrix(0,nrow=(nrow(data)-1),ncol=(ncol(data)-1))
    for (curstep in c(1:numsteps))
    {
        if (curstep == numsteps)
        {
            for (curcol in c(2:ncol(data)))
            {
                # evaluate \Gamma vector across all samples at a particular time
                # this is a matrix because we're also evaluating at xvec
                gammamat = integrandmat(data[2:nrow(data),curcol],xvec,h,driftfun,difffun,c0)
                likelihood[,(curcol-1)] = k*gammamat %*% pdfmatrix[,(curcol-1)]

                gd1mat = Dtheta1(data[2:nrow(data),curcol],xvec,h,driftfun,difffun1,c0)
                gd2mat = Dtheta2(data[2:nrow(data),curcol],xvec,h,driftfun,difffun1,c0)
                d1mat[,(curcol-1)] = k*gd1mat %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta1[,(curcol-1)]
                d2mat[,(curcol-1)] = k*gd2mat %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta2[,(curcol-1)]
            }
        }
        pdfmatrix = k*A %*% pdfmatrix
        qmattheta1 = k*A %*% qmattheta1 + k*D1 %*% pdfmatrix
        qmattheta2 = k*A %*% qmattheta2 + k*D2 %*% pdfmatrix
    }

    return(list(x=xvec,y=pdfmatrix,z1=qmattheta1,z2=qmattheta2,lik=likelihood,likd1=d1mat,likd2=d2mat))
}


