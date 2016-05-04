rm(list=ls(all=TRUE))

library('nloptr')
load('fakedata.RData')
fd = xtraj
source('adjoint.R')

objgradfun <- function(theta)
{
	out = gradobjfun(c0=theta, h=myh, k=myk, bigm=mybigm, littlet=littletnew, data=fd)
	return(list("objective" = out$obj, "gradient" = out$grad))
}

myh = 0.05
myk = (myh)^(0.75)
mybigm = ceiling(pi/(myk^1.5))
theta = {}
theta[1:4] = c(0, 0, 0, 0)
theta[5:7] = c(0, 0, 0)
theta[8] = 0.7
littletnew=0.5


res <- nloptr(x0 = theta, eval_f = objgradfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))


