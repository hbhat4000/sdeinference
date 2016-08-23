rm(list = ls(all = TRUE))

# load data
load('fakedata.RData')
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
        gradient[i] = -sum(probmat$grad[[i]] / probmat$lik)

    return(list("objective"=objective,"gradient"=gradient))
}

# create grid, compute densities on that grid
myh = 0.001
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.1))

# initial condition fakedata = c(1,4,0.5)
theta = c(rep(0,4),0.5)

library('nloptr')

res <- nloptr(x0 = theta, eval_f = objgradfun, lb = c(rep(-10,4),0.2), ub = c(rep(10,4),1), opts = list("algorithm"="NLOPT_LD_MMA", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))

save(res,file='optimtestresult.RData')


