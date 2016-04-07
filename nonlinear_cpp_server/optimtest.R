rm(list = ls(all = TRUE))

# library for optimization algorithms
library('nloptr')

# load the dtq package for finding the computing likelihood function
library('Rgdtq')

# load data
load('fakedata.RData')
fd = xtraj[2:nrow(xtraj),]

objgradfun <- function(theta)
{
    numparam = length(theta)
    # OLD: probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)
    # probmat = array(0, dim = c(nrow(fd), ncol(fd) - 1, numparam + 1))
    probmat = Rgdtq(thetavec = theta, h = myh, k = myk, M = mybigm, littlet = 1, init_data = fd)

    mylik = probmat[,,1]

    mylik[mylik < 0] = 0

    objective = -sum(log(mylik))
    gradient = numeric(length = numparam)

    for (i in c(1:numparam))
        gradient[i] = -sum(probmat[,,i] / mylik)

    return(list("objective" = objective, "gradient" = gradient))
}

# create grid, compute densities on that grid
myh = 0.05
myk = (myh)^(0.75)
mybigm = ceiling(pi/(myk^1.5))

theta = c(2, 2, 0.7)

res <- nloptr(x0 = theta, eval_f = objgradfun, lb = c(0.1, 0, 0.1), ub = c(10, 10, 1.0), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-3))

