rm(list = ls(all = TRUE))

# load data
load('fakedata.RData')
fd = xtraj[1:4,1:4]

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
myh = 0.1
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.5))

# initial condition fakedata = c(1,4,0.5)
theta = c(1, 0.5, 0.5)

#library('nloptr')

#res <- nloptr(x0 = theta, eval_f = objgradfun, lb = c(0.1, 0, 0.1), ub = c(4, 4, 4), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))

#fname = paste("optimtest2_",myh,".RData",sep='')
#save(res,file=fname)

