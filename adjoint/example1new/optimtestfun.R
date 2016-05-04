library('nloptr')
load('fakedata.RData')
tvec = xtraj[1,]
fd = xtraj[-1,]
source('adjoint.R')

myh = 0.1
myk = (myh)^(0.6)
mybigm = ceiling(pi/(myk^1.5))
numparam = 6
inittheta = as.numeric(0*c(1:(numparam-1)))
inittheta[numparam] = 0.7

myoptfun <- function(data, indices)
{
	fdnewdata = rbind(tvec, fd[indices, ])
	objgradfun <- function(theta)
	{
		out = gradobjfun(c0=theta, h=myh, k=myk, bigm=mybigm, littlet=1, data=fdnewdata)
		return(list("objective" = out$obj, "gradient" = out$grad))
	}


	res <- nloptr(x0 = inittheta, eval_f = objgradfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))
	return(res)
}
