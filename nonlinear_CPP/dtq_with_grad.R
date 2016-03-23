driftfun <- function(c0, y) 
{
  # computes drift as a matrix
  f = (c0[1])*(y)*(c0[2] - y^2)
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

cdt <- function(c0, h, k, bigm, littlet, pdfmatrix0, data)
{
    numsteps = ceiling(littlet/h)
    xvec = c((-bigm):bigm)*k
    npts = length(xvec)

    A = integrandmat(xvec, xvec, h, driftfun, difffun, c0)
    D = Dtheta(xvec, xvec, h, driftfun, difffun1, c0)
    
    # derivatives of the pdf w.r.t theta_vec
    # qmattheta = list(0*pdfmatrix0, 0*pdfmatrix0, 0*pdfmatrix0)
    qmattheta1 = 0*pdfmatrix0
    qmattheta2 = 0*pdfmatrix0
    qmattheta3 = 0*pdfmatrix0
    
    # evolve all pdf's forward by correct # of time steps
    pdfmatrix = pdfmatrix0
    likelihood = matrix(0, nrow = (nrow(data) - 1), ncol = (ncol(data)-1))
    #dmat = list(0*likelihood, 0*likelihood, 0*likelihood)
    d1mat = 0*likelihood
    d2mat = 0*likelihood
    d3mat = 0*likelihood
    
    for (curstep in c(1:numsteps))
    {
        if (curstep == numsteps)
        {
            for (curcol in c(2:ncol(data)))
            {
                # evaluate \Gamma vector across all samples at a particular time
                # this is a matrix because we're also evaluating at xvec
              
                gammamat = integrandmat(data[2:nrow(data), curcol], xvec, h, driftfun, difffun, c0)
                likelihood[,(curcol-1)] = k*gammamat %*% pdfmatrix[,(curcol-1)]

                gdmat = Dtheta(data[2:nrow(data),curcol], xvec, h, driftfun, difffun1, c0)

                d1mat[,(curcol-1)] = k*gdmat[[1]] %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta1[,(curcol-1)]
                d2mat[,(curcol-1)] = k*gdmat[[2]] %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta2[,(curcol-1)]
                d3mat[,(curcol-1)] = k*gdmat[[3]] %*% pdfmatrix[,(curcol-1)] + k*gammamat %*% qmattheta3[,(curcol-1)]
            }
        }

        pdfmatrix = k*A %*% pdfmatrix
        qmattheta1 = k*A %*% qmattheta1 + k*D[[1]] %*% pdfmatrix
        qmattheta2 = k*A %*% qmattheta2 + k*D[[2]] %*% pdfmatrix
        qmattheta3 = k*A %*% qmattheta3 + k*D[[3]] %*% pdfmatrix
    }

    return(list(x = xvec, y = pdfmatrix, z1 = qmattheta1, z2 = qmattheta2, z3 = qmattheta3, lik = likelihood, likd1 = d1mat, likd2 = d2mat, likd3 = d3mat))
}


