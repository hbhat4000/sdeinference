rm(list = ls(all = TRUE))

# library for optimization algorithms
library('nloptr')

# load the dtq package for finding the computing likelihood function
library('Rgdtq')

# load data
load('fakedata.RData')
fd = xtraj[[2]]

# load necessary functions
# source('dtq_with_grad.R')

objgradfun <- function(theta)
{
    # OLD: probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)
    probmat = Rgdtq(thetavec = theta, h = myh, k = myk, M = mybigm, littlet = 1, init_data = fd)
    
    nrowmat = nrow(probmat) / 4
    ncolmat = ncol(probmat)
    mylik = probmat[1:nrowmat,1:ncolmat]
    grad1 = probmat[(nrowmat+1):(2*nrowmat),1:ncolmat]
    grad2 = probmat[(2*nrowmat+1):(3*nrowmat),1:ncolmat]
    grad3 = probmat[(3*nrowmat+1):(4*nrowmat),1:ncolmat]

    numparam = length(theta)

    # mylik = probmat$lik
    mylik[mylik < 0] = 0

    objective = -sum(log(mylik))
    gradient = numeric(length = numparam)

    gradient[1] = -sum(grad1 / mylik)
    gradient[2] = -sum(grad2 / mylik)
    gradient[3] = -sum(grad3 / mylik)

    # for (i in c(1:nc0))
    #     gradient[i] = -sum(probmat$grad[[i]] / probmat$lik)

    return(list("objective" = objective, "gradient" = gradient))
}

# create grid, compute densities on that grid
myh = 0.01
myk = (myh)
mybigm = ceiling(pi/(myk^1.5))

theta = c(2, 2, 0.7)

res <- nloptr(x0 = theta, eval_f = objgradfun, lb = c(0.1, 0, 0.1), ub = c(10, 10, 1.0), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-3))

