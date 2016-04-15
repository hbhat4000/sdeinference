rm(list = ls(all = TRUE))

# library for optimization algorithms
library('nloptr')
library('Rcpp')
library('RcppArmadillo')

library('Rgdtq')

# load the dtq file for computing the likelihood function
# source("rcpp_dtq.cpp")

# load data
load('fakedata.RData')
fd = xtraj

f <- function(z)
{
    theta = z$thetavec
    y = as.matrix(z$Y)
    # return(theta[1] * (theta[2] - y))
    out = matrix(5,nrow=3,ncol=3)
    return(out)
}

g <- function(z)
{
    theta = z$thetavec
    out = matrix(theta[3],nrow=nrow(z$Y),ncol=ncol(z$Y))
    # return(out)
    out = matrix(5,nrow=3,ncol=3)
    return(out)
}


objgradfun <- function(theta)
{
    numparam = length(theta)
    # probmat = array(0, dim = c(nrow(fd), ncol(fd) - 1, numparam + 1))
    probmat = gDTQ(thetavec = theta, h = myh, k = myk, M = mybigm, littlet = 1, init_data = fd, drift=f, diffusion=g)

    mylik = probmat[,,1]
    mylik[mylik < 0] = 0

    objective = -sum(log(mylik))
    gradient = numeric(length = numparam)

    for (i in c(1:numparam))
        gradient[i] = -sum(probmat[,,(i+1)] / mylik)

    return(list("objective" = objective, "gradient" = gradient))
}

# create grid, compute densities on that grid
myh = 0.1
myk = (myh)^(0.75)
mybigm = ceiling(pi/(myk^1.5))

theta = c(2, 1, 1)

res <- nloptr(x0 = theta, eval_f = objgradfun, lb = c(0.1, 0, 0.1), ub = c(10, 10, 1.0), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-3))

