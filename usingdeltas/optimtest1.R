rm(list=ls(all=TRUE))


# load data
load('fakedata_ou_noise.RData')
fd = xtraj

# load necessary functions
source('cdt_new.R')

# function that computes log likelihood
mylikelihood <- function(likden,dat)
{
      numsteps = dim(likden$y)[2]
	logpost = 0
	for (k in c(1:(numsteps-1)))
	{
		# browser()
    		likdat = approx(x=likden$x,y=likden$y[,k],xout=dat[2,k+1])
    		logpost = log(likdat$y) + logpost	
	}
return(logpost)
}

myh = 0.1

objectfun <- function(c0)
{
	probmat = cdt(c0,h=myh,littlet=1,data=fd)
	-mylikelihood(likden=probmat,dat=fd)
}
optim(c(0.8, 0.9), objectfun)