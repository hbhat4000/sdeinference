rm(list=ls(all=TRUE))

library('nloptr')
load('fakedata_t1_T50_nt100.RData')
fd = xtraj
dimvec = dim(fd)
vec1 = c(1:ceiling(dimvec[2]/2))
fdtrain=fd[,vec1]
fdtest = fd[-1,-vec1]

source('adjoint.R')

objgradfun <- function(theta)
{
	out = gradobjfun(c0=theta, h=myh, k=myk, bigm=mybigm, littlet=1, data=fdtrain,  diffcoeff)
	return(list("objective" = out$obj, "gradient" = out$grad))
}

myh = 0.1
myk = (myh)^(0.6)
mybigm = ceiling(pi/(myk^1.5))
theta = as.numeric(0*c(1:5))
diffcoeff=1

res <- nloptr(x0 = theta, eval_f = objgradfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))

inverted_theta = res$solution
print(inverted_theta)
source('driftcomparison.R')
mse_path = compare (thetaapprox=inverted_theta, fdtest, fdtrain, diffcoeff)