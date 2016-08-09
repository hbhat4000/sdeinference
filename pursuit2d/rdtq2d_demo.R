# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data which has 2 lists, chaser and runner, each with t, x, y
load('fakedata.RData')
tchase = chaser[[1]]
xchase = chaser[[2]]
ychase = chaser[[3]]
trun = runner[[1]]
xrun = runner[[2]]
yrun = runner[[3]]

truethetavec = c(1.2, 1, 1)

# algorithm parameters
mydatapoints = length(tchase) - 1
myh = 0.05
myk = (0.8*myh)^0.75
xylimit = 94    # court dimension is 94*50
# xylimit = max(abs(xchase), abs(ychase))

# check PDF
# mycheck = PDFcheck(thetavec = truethetavec, h = myh, k = myk, yM = xylimit)
# print(c(min(mycheck), mean(mycheck), max(mycheck)))

xchase = xchase[1:mydatapoints]
ychase = ychase[1:mydatapoints]
xrun = xrun[1:mydatapoints]
yrun = yrun[1:mydatapoints]

# Metropolis Hastings parameters
burnin = 100
numsteps = 1000
totsteps = numsteps + burnin
mcmc = numeric(length = totsteps)
mcmc[1] = 0.1

# time increment from data
timeinc = tchase[2] - tchase[1]
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
    datamat = matrix(c(dat[[2]], dat[[3]]), ncol = 2)
    logpost = 0

    for(i in c(1:steps))
    { 
      myloc = datamat[(i+1),c(1:2)]
      likdat = interp.surface(obj = list(x = xvec, y = xvec, z = matrix(likden[,i], nrow = mm)), loc = myloc)
      likdat[likdat <= 2.2e-16] = 2.2e-16
      logpost = logpost + sum(log(likdat))
    }
    logpost = logpost + prior
    return(logpost)
}

thetavec = truethetavec
thetavec[1] = mcmc[1]
oldden = Rdtq2d(thetavec, xchase, ychase, xrun, yrun, h = myh, numsteps = myns, k = myk, yM = xylimit)

# checkk = PDFcheck(thetavec, h = myh, k = myk, yM = xylimit)
# print(checkk)

oldpost = myposterior(likden = oldden, dat = chaser, prior = myprior(mcmc[1]))
artrack = numeric(length = (totsteps-1))

for (i in c(1:(totsteps-1)))
{
    # generate proposal
    z = rnorm(n = 1, sd = 0.25)
    prop = mcmc[i] + z
	
    # calculate likelihood and posterior density
    thetavec = truethetavec
    thetavec[1] = prop
    propden = Rdtq2d(thetavec, C1, C2, h = myh, numsteps = myns, k = myk, yM = xylimit)
    proppost = myposterior(likden = propden, dat = chaser, prior = myprior(prop)) 
    rho = exp(proppost - oldpost)
    maxcolsumerr = max(abs(colSums(propden)*myk^2 - 1))

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        mcmc[i+1] = prop
        oldden = propden
        oldpost = proppost
        print(paste("Accepted the proposal", prop))
        artrack[i] = 1
    }
    else
    {
        mcmc[i+1] = mcmc[i]
        print(paste("Rejected the proposal", prop))
        artrack[i] = 0
    }
    print(c(i, rho, maxcolsumerr, mcmc[i+1]))
    flush.console()
}

# throw away the burnin steps
mcmc = mcmc[(burnin+1):totsteps]

# save everything
save.image(file = 'posterior1.RData')

# keep every 10th step
mcmc = mcmc[seq(from = 1, to = numsteps, by = 10)]

# plot kernel density estimate of samples
# overlay true posterior in red
myden = density(mcmc, bw = 'SJ')
plot(myden)

