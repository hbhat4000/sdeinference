rm(list=ls(all=TRUE))

library('nloptr')
load('fakedata.RData')
fd = xtraj
source('adjointpar.R')

objgradfun <- function(theta)
{
	out = gradobjfun(c0=theta, h=myh, k=myk, bigm=mybigm, littlet=1, data=fd,  diffcoeff=1)
	return(list("objective" = out$obj, "gradient" = out$grad))
}

myh = 0.02
myk = (myh)^(0.6)
mybigm = ceiling(pi/(myk^1.5))
numpara = 6
theta = as.numeric(0*c(1:numpara))
theta[numpara] = 0.7


myfunc <- function()
{
    #mybounds = numeric(length=length(theta))
    #mybounds[seq(from=1,to=length(theta),by=2)] = 0.1
    #mybounds[seq(from=2,to=length(theta),by=2)] = 25
    #res <- nloptr(x0 = theta, eval_f = objgradfun, lb=-mybounds, ub=mybounds, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))
    res <- nloptr(x0 = theta, eval_f = objgradfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))
    return(res)
}

# to recover Hessian
result = myfunc()
theta = result$solution

objfun <- function(theta) { return(objgradfun(theta)$objective) }
gradfun <- function(theta) { return(objgradfun(theta)$gradient) }
invhess = solve(optimHess(par=theta,fn=objfun,gradfun))
diagentries = diag(invhess)
diagentries[diagentries < 0] = 0
print(sqrt(diagentries))



