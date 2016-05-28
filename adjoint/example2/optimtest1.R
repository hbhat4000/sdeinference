rm(list=ls(all=TRUE))

library('nloptr')
load('fakedata.RData')
fd = xtraj
source('adjoint.R')

objgradfun <- function(theta)
{
	out = gradobjfun(c0=theta, h=myh, k=myk, bigm=mybigm, littlet=1, data=fd,  diffcoeff=1)
	return(list("objective" = out$obj, "gradient" = out$grad))
}

myh = 0.1
myk = (myh)^(0.6)
mybigm = ceiling(pi/(myk^1.5))
# theta = as.numeric(c(1:8))
theta = rnorm(12)

res <- nloptr(x0 = theta, eval_f = objgradfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))


