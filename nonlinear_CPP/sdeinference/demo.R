rm(list=ls(all=TRUE))

# load data
load('fakedata_ou_noise.RData')
fd = xtraj

# function that computes log likelihood
mylikelihood <- function(likden,dat)
{
    numsteps = dim(likden$y)[2]
	logpost = 0
	for (k in c(1:(numsteps-1)))
	{
    	likdat = approx(x=likden$x,y=likden$y[,k],xout=dat[2,k+1])
    	logpost = log(likdat$y) + logpost	
	}
return(logpost)
}

gk <- function(likden, dat)
{
     numsteps = dim(likden$z1)[2]
	val1 = 0
	val2 = 0
	for (k in c(1:(numsteps-1)))
	{
		
    	prob = approx(x=likden$x,y=likden$y[,k],xout=dat[2,k+1])
		diffprob1 = approx(x=likden$x,y=likden$z1[,k],xout=dat[2,k+1])
		diffprob2 = approx(x=likden$x,y=likden$z2[,k],xout=dat[2,k+1])
    	val1 = diffprob1$y/prob$y + val1
		val2 = diffprob2$y/prob$y + val2
	}
	vec = c(val1, val2)
return(vec)
}

myh = 0.1

objectfun <- function(c0)
{
	probmat = cdt(c0,h=myh,littlet=1,data=fd)
	-mylikelihood(likden=probmat,dat=fd)
}

gradientfun <- function(c0)
{
	probmat = cdt(c0,h=myh,littlet=1,data=fd)
	-gk(likden = probmat, dat = fd)
}

# true parameters
theta = c(0.9, 1)

optim(c(0.8,0.9), objectfun)
# [1] 0.9298092 0.9526505



(res <- optim(c(0.8,0.9), objectfun, gradientfun, method = "BFGS"))
# [1] 0.8616467 0.9616467
optimHess(res$par, objectfun, gradientfun)
#         [,1]      [,2]
#[1,] 66.30338 31.723790
#[2,] 31.72379 -2.855804


# optim(c(0.8,0.9), objectfun,  gradientfun, method = "CG")
#  0.8617357 0.9617357


#optim(c(0.8,0.9), objectfun,  gradientfun, method = "BFGS")
# 0.8616467 0.9616467