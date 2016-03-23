# clear all memory
rm(list=ls(all=TRUE))

# load the dtq package for finding the optimized likelihood function
library('Rgdtq')

# load package for interpolation
# library('fields')

# load fake data in X which has x1, x2, t
load('fakedata.RData')

# load true thetavec
source('truethetavec.R')

# algorithm parameters
mydatapoints = nrow(X) - 1
myh = 0.005
powk = 0.5
myk = 0.8*(myh^powk)
xylimit = 3*max(abs(X[,1]))
C0 = X[1:mydatapoints,1]
burnin = 100
numsteps = 1000
totsteps = numsteps + burnin
infer = numeric(length = totsteps)
infer[1] = 1.0

thetavec = truethetavec
thetavec[1] = infer[1]

# checkk = PDFcheck(thetavec,h=myh,k=myk,yM=xylimit)
# print(sum(checkk > 0.9))
# print(length(checkk))
# acceptablenorm = sum(checkk > 0.9)/(length(checkk))
# print(c(powk, myk, acceptablenorm))
# flush.console()
# while(acceptablenorm < 0.95)
# #while(min(checkk)<0.8)
# {
#     powk = powk*(1.1)
#     myk = myh^powk
#     checkk = PDFcheck(thetavec,h=myh,k=myk,yM=xylimit)
#     print(sum(checkk > 0.9))
#     print(length(checkk))
#     acceptablenorm = sum(checkk > 0.9)/(length(checkk))
#     print(c(powk, myk, acceptablenorm))
#     flush.console()
# }

# time increment from data
timeinc = X[2,2] - X[1,2]
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
# mm = length(xvec)

myposterior <- function(likden, dat, prior)
{
    # for scenario 2, likden has 'datapoints' pdfs  
    logpost = 0
    for(i in c(1:steps))
    { 
      likdat = approx(x = xvec, y = dat[i,1])
      likdat[likdat <= 2.2e-16] = 2.2e-16
      logpost = logpost + sum(log(likdat))
    }

    logpost = logpost + prior
    return(logpost)
}

oldden = Rgdtq(thetavec, C0, h = myh, numsteps = myns, k = myk, yM = xylimit)
oldpost = myposterior(likden = oldden, dat = X, prior = myprior(infer[1]))
artrack = numeric(length = (totsteps-1))

for (i in c(1:(totsteps-1)))
{
    # generate proposal
    z = rnorm(n = 1, sd = 0.25)
    prop = infer[i] + z
	
    # calculate likelihood and posterior density
    thetavec = truethetavec
    thetavec[1] = prop
    propden = Rgdtq(thetavec, C0, h = myh, numsteps = myns, k = myk, yM = xylimit)
    proppost = myposterior(likden = propden, dat = X, prior = myprior(prop)) 
    rho = exp(proppost - oldpost)
    maxcolsumerr = max(abs(colSums(propden)*myk^2 - 1))

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        infer[i+1] = prop
        oldden = propden
        oldpost = proppost
        print(paste("Accepted the proposal",prop))
        artrack[i] = 1
    }
    else
    {
        infer[i+1] = infer[i]
        print(paste("Rejected the proposal",prop))
        artrack[i] = 0
    }
    print(c(i, rho, maxcolsumerr, infer[i+1]))
    flush.console()
}

# throw away the burnin steps
infer = infer[(burnin + 1):totsteps]

# save everything
save.image(file = 'posteriorsamples.RData')

# keep every 10th step
infer = infer[seq(from = 1, to = numsteps, by = 10)]

# plot kernel density estimate of samples
# overlay true posterior in red
myden = density(infer,bw='SJ')
plot(myden)

