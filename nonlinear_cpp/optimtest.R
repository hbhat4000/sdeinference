rm(list = ls(all = TRUE))

# library for optimization algorithms
library('nloptr')

# load the dtq package for finding the computing likelihood function
library('Rgdtq')

# load required functions and parameters
source('parameters.R')

# load data which is saved as fd 
fname = paste('fakedata_', fakedatanum, '.RData', sep = '')
load(fname)
fd = xtraj

objgradfun <- function(theta)
{
    numparam = length(theta)
    # probmat = array(0, dim = c(nrow(fd), ncol(fd) - 1, numparam + 1))
    probmat = Rgdtq(thetavec = theta, h = gridh, k = gridk, M = gridM, littlet = 1, init_data = fd)

    mylik = probmat[,,1]
    mylik[mylik < 0] = 0

    objective = -sum(log(mylik))
    gradient = numeric(length = numparam)

    for (i in c(1:numparam))
        gradient[i] = -sum(probmat[,,(i+1)] / mylik)

    return(list("objective" = objective, "gradient" = gradient))
}

print(objgradfun(theta))
# res <- nloptr(x0 = theta, eval_f = objgradfun, lb = c(0.1, 0, 0.1), ub = c(10, 10, 2), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4))