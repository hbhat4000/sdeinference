rm(list=ls(all=TRUE))


# load data
load('fakedata_ou_noise.RData')
fd = xtraj

# load necessary functions
source('qdt_with_grad.R')

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

gkfun <- function(likden, dat)
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
	-gkfun(likden = probmat, dat = fd)
}


c0 = c(0.5,0.7)
fk = objectfun(c0)
gk = gradientfun(c0)
normerror = 1
alpha = 1
c = 10^(-4)

#steepest descent

while(normerror > 0.01)
{
	newc0 = c0-alpha*gk
      fknew = objectfun(newc0)
	
 	if ( fknew < fk + c*alpha*(t(gk) %*% gk) )
	{
		c0 = newc0
		gk = gradientfun(c0)
		normerror = sqrt(sum(gk^2))
		fk = fknew	
	}else
	{
 	alpha = 0.5*alpha 
	}	
}


print(normerror)
print(c0)