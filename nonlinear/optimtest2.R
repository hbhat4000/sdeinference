rm(list = ls(all = TRUE))

library('KernSmooth')

# load data
load('fakedata_ou_noise.RData')
fd = xtraj

# load necessary functions
source('dtq_with_grad.R')

objectfun <- function(c0)
{
    probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, pdfmatrix = ics, data = fd)
    mylik = probmat$lik
    mylik[mylik < 0] = 0
    out = -sum(log(mylik))
    return(out)
}

gradientfun <- function(c0)
{
	probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, pdfmatrix = ics, data = fd)
    out = c(0,0,0)
        
    out[1] = -sum(probmat$likd1 / probmat$lik)
    out[2] = -sum(probmat$likd2 / probmat$lik)
    out[3] = -sum(probmat$likd3 / probmat$lik)
    return(out)
}

# create grid, compute densities on that grid
myh = 0.05
myk = myh^0.95
mybigm = ceiling(pi/(myk^1.5))
nsp = nrow(fd) - 1
nts = ncol(fd) - 1
ics = matrix(0, nrow = (2*mybigm + 1), ncol = nts)
for (i in c(1:nts))
{
    myden = bkde(x = fd[2:(nsp+1),i], gridsize = (2*mybigm + 1), range.x = c(-mybigm*myk,mybigm*myk))
    ics[,i] = myden$y
}

# initial condition fakedata = c(1,4,0.5)
theta = c(1, 7, 0.5)

library('nloptr')

res <- nloptr(x0 = theta, eval_f = objectfun, eval_grad_f = gradientfun, lb = c(0, 0, 0.1), ub = c(10, 10, 1.0), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_rel"=1e-8))
#res <- nloptr(x0 = theta, eval_f = objectfun, eval_grad_f = gradientfun, lb = c(0, 0, 0.1), ub = c(10, 10, 1.0), opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = TRUE, "xtol_rel"=1e-8))
