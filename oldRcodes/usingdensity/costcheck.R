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
    		logpost = log(likdat$y) +logpost
    		
	}
return(logpost)
}


c1vec = seq(from = 0.4, to=1.4, by=0.1)
c2vec = seq(from = 0.2, to=1.5, by=0.1)

myh = 0.05

n1 = length(c1vec)
n2 = length(c2vec)

cost = matrix(0, nrow = n1, ncol = n2)
for (i in c(1:n1))
{
	for (j in c(1:n2))
	{
	print(c(i,j))
	c0 = c(c1vec[i], c2vec[j])
	probmat = cdt(c0,h=myh,littlet=1,data=fd)
	cost[i,j] = -mylikelihood(likden=probmat,dat=fd)
	}
}
minval = min(cost)-10
maxval = max(cost)
if (maxval == Inf)
{
	maxval = 3000
}
contour(c1vec,c2vec,cost,seq(floor(minval),ceiling(maxval) , 10))


print(min(cost))

indx1 = which(c1vec==0.9)
indx2 = which(c2vec==1)
print(cost[indx1,indx2])
