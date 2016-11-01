# clear all memory
rm(list=ls(all=TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data in X which has x1, x2, t
load('fakedata_nonlinear.RData')

# keep every 100th row
X = X[seq(from=1,to=nrow(X),by=100),]

# load true thetavec
source('truethetavec.R')

# algorithm parameters
# time increment from data and time step
mydatapoints = nrow(X) - 1
timeinc = X[2,3] - X[1,3]
myh = timeinc/4
myns = floor(timeinc/myh)
print(c(myh,myns))
myk = 0.8*myh^0.75
xylimit = 2*max(abs(X[,1:2]))

C1 = X[1:mydatapoints,1]
C2 = X[1:mydatapoints,2]
burnin = 100
numsteps = 1000
totsteps = numsteps + burnin

x = matrix(0, nrow = totsteps, ncol = 2)
x[1,] = c(0.1,0.1)

# define log prior
myprior <- function(z)
{
    return(dnorm(x = z, mean = 0, sd = 100, log = TRUE))
}

# function that computes log posterior

M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)
mm = length(xvec)

myposterior <- function(likden, dat, prior)
{
    steps = ncol(likden)
    # for scenario 2, likden has 'datapoints' pdfs  
    logpost = 0
    for(i in c(1:steps))
    { 
      myloc = matrix(dat[(i+1),c(1:2)],ncol=2)
      likdat = interp.surface(obj=list(x=xvec,y=xvec,z=matrix(likden[,i],nrow=mm)),loc=myloc)
      likdat[likdat <= 2.2e-16] = 2.2e-16
      logpost = logpost + sum(log(likdat))
    }
    logpost = logpost + sum(prior)
    return(logpost)
}

thetavec = truethetavec
thetavec[1] = x[1,1]
thetavec[2] = x[1,2]
oldden = Rdtq2d(thetavec,C1,C2,h=myh,numsteps=myns,k=myk,yM=xylimit)
checkk = PDFcheck(thetavec,h=myh,k=myk,yM=xylimit)
print(c(min(checkk),mean(checkk),max(checkk)))
oldpost = myposterior(likden=oldden, dat=X, prior=myprior(x[1,]))
artrack = numeric(length=(totsteps-1))

for (i in c(1:(totsteps-1)))
{
    # generate proposal
    z = rnorm(n = 2, sd = 0.1)
    prop = x[i,] + z
	
    # calculate likelihood and posterior density
    thetavec = truethetavec
    thetavec[1] = prop[1]
    thetavec[2] = prop[2]
    propden = Rdtq2d(thetavec,C1,C2,h=myh,numsteps=myns,k=myk,yM=xylimit)
    proppost = myposterior(likden=propden, dat=X, prior=myprior(prop)) 
    rho = exp(proppost-oldpost)
    maxcolsumerr = max(abs(colSums(propden)*myk^2 - 1))

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        x[i+1,] = prop
        oldden = propden
        oldpost = proppost
        print(paste("Accepted the proposal",prop))
        artrack[i] = 1
    }
    else
    {
        x[i+1,] = x[i,]
        print(paste("Rejected the proposal",prop))
        artrack[i] = 0
    }
    print(c(i,rho,maxcolsumerr,x[i+1,]))
    flush.console()
}

# throw away the burnin steps
x = x[(burnin+1):totsteps,]

# save everything
save.image(file = 'ps_both.RData')

# percentage difference between empirical and true means
print((mean(x[,1]^2) - 2*pi)/(2*pi))

# percentage difference between empirical and true modes
myden = density(x[,1]^2)
print((myden$x[which.max(myden$y)] - 2*pi)/(2*pi))


