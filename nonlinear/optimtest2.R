rm(list = ls(all = TRUE))

# load data
load('fakedata5.RData')
fd = xtraj

# load necessary functions
source('dtq_with_grad.R')

objgradfun <- function(c0)
{
    probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)
    mylik = probmat$lik
    mylik[mylik < 0] = 0
    objective = -sum(log(mylik))
    nc0 = length(c0)
    gradient = numeric(length=nc0)
    for (i in c(1:nc0))
        gradient[i] = - sum(probmat$grad[[i]] / probmat$lik)

    return(list("objective"=objective,"gradient"=gradient))
}

# create grid, compute densities on that grid
myh = 0.01
myk = myh
mybigm = ceiling(pi/(myk^1.5))

# initial condition fakedata = c(1,4,0.5)
theta = c(4, 4, 0.5)

library('nloptr')

res <- nloptr(x0 = theta, eval_f = objgradfun, lb = c(0, 0, 0.01), ub = c(10, 10, 2), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_abs"=1e-4))

