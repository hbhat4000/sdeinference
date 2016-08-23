# clear all memory
rm(list=ls(all=TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data in X which has x1, x2, t
load('fakedata1.RData')

# delete every other row
X = X[seq(from=1,to=nrow(X),by=2),]

# load true thetavec
source('truethetavec.R')
truethetavec[2] = sqrt(truethetavec[2])

# algorithm parameters
mydatapoints = nrow(X) - 1
myh = 0.005
myk = 0.8*myh^0.75
xylimit = 2*max(abs(X[,1:2]))

# check PDF
mycheck = PDFcheck(thetavec=truethetavec,h=myh,k=myk,yM=xylimit)
print(c(min(mycheck),mean(mycheck),max(mycheck)))

C1 = X[1:mydatapoints,1]
C2 = X[1:mydatapoints,2]
burnin = 100
numsteps = 1000
totsteps = numsteps + burnin
x = numeric(length = totsteps)
x[1] = 0.1

# time increment from data
timeinc = X[2,3] - X[1,3]
myns = floor(timeinc/myh)
print(myns)

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
    logpost = logpost + prior
    return(logpost)
}

thetavec = truethetavec
thetavec[1] = x[1]
oldden = Rdtq2d(thetavec,C1,C2,h=myh,numsteps=myns,k=myk,yM=xylimit)
checkk = PDFcheck(thetavec,h=myh,k=myk,yM=xylimit)
print(checkk)
oldpost = myposterior(likden=oldden, dat=X, prior=myprior(x[1]))
artrack = numeric(length=(totsteps-1))

for (i in c(1:(totsteps-1)))
{
    # generate proposal
    z = rnorm(n = 1, sd = 0.25)
    prop = x[i] + z
	
    # calculate likelihood and posterior density
    thetavec = truethetavec
    thetavec[1] = prop
    propden = Rdtq2d(thetavec,C1,C2,h=myh,numsteps=myns,k=myk,yM=xylimit)
    proppost = myposterior(likden=propden, dat=X, prior=myprior(prop)) 
    rho = exp(proppost-oldpost)
    maxcolsumerr = max(abs(colSums(propden)*myk^2 - 1))

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        x[i+1] = prop
        oldden = propden
        oldpost = proppost
        print(paste("Accepted the proposal",prop))
        artrack[i] = 1
    }
    else
    {
        x[i+1] = x[i]
        print(paste("Rejected the proposal",prop))
        artrack[i] = 0
    }
    print(c(i,rho,maxcolsumerr,x[i+1]))
    flush.console()
}

# throw away the burnin steps
x = x[(burnin+1):totsteps]

# save everything
save.image(file = 'posteriorsamples_val3.RData')

# keep every 10th step
x = x[seq(from = 1, to = numsteps, by = 10)]

# plot kernel density estimate of samples
# overlay true posterior in red
myden = density(x,bw='SJ')
plot(myden)

