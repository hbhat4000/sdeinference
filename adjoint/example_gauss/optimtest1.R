rm(list=ls(all=TRUE))

library('nloptr')
load('fakedata_pointfive_100.RData')
fd = xtraj
source('adjoint.R')

objgradfun <- function(theta)
{
	out = gradobjfun(c0=theta, h=myh, k=myk, bigm=mybigm, littlet=littletnew, data=fd)
	return(list("objective" = out$obj, "gradient" = out$grad))
}

myh = 0.005
myk = (myh)^(0.6)
mybigm = ceiling(pi/(myk^1.5))
npara = 12
theta = as.numeric(0*c(1:(npara-1)))
theta[npara] = 0.7
littletnew = 0.5

res <- nloptr(x0 = theta, eval_f = objgradfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))


