# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data which has 2 lists, chaser and runner, each with t, x, y
load('fakedata3.RData')
tchase = chaser[[1]]
xchase = chaser[[2]]
ychase = chaser[[3]]
runner = matrix(unlist(runner),ncol=3)

truethetavec = c(0.3, 0.5)

# algorithm parameters
mydatapoints = length(tchase) - 1

myh = 0.4
# myk = 0.8*(myh)^0.9
myk = 0.5
xylimit = 5    # court dimension is 94*50

# xylimit = max(abs(xchase), abs(ychase))

# check PDF
# mycheck = PDFcheck(thetavec = truethetavec, h = myh, k = myk, yM = xylimit)
# print(c(min(mycheck), mean(mycheck), max(mycheck)))

xchase = xchase[1:mydatapoints]
ychase = ychase[1:mydatapoints]

# Metropolis Hastings parameters
burnin = 10
numsteps = 100
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
    # print(steps)
    datamat = matrix(c(dat[[2]], dat[[3]]), ncol = 2)
    # print(datamat)
    logpost = 0

    for(i in c(1:steps))
    { 
      myloc = matrix(datamat[(i+1),c(1:2)], ncol = 2)
      # print(myloc)

      # interp.surface(obj, loc)
      # obj : a list with components x, y and z s.t. x and y are the X and Y grid values and z is a matrix
      # with the corresponding values of the surface
      # loc : a matrix of irregular locations to interpolate, first column of loc is the X coordinates
      # and second column is the Y's
      likdat = interp.surface(obj = list(x = xvec, y = xvec, z = matrix(likden[,i], nrow = mm)), loc = myloc)
      likdat[likdat <= 2.2e-16] = 2.2e-16
      logpost = logpost + sum(log(likdat))
    }
    logpost = logpost + prior
    return(logpost)
}

thetavec = truethetavec
thetavec[1] = mcmc[1]

# print("Pass 1")

# gammavec = rep(1,mydatapoints+1)
gammavec = c(1, 0.8)

oldden = Rdtq2d(thetavec, gammavec, runner, c1=xchase, c2=ychase, h = myh, numsteps = myns, k = myk, yM = xylimit)

# print("Pass 2")
# checkk = PDFcheck(thetavec, h = myh, k = myk, yM = xylimit)
# print(checkk)

oldpost = myposterior(likden = oldden, dat = chaser, prior = myprior(mcmc[1]))
artrack = numeric(length = (totsteps-1))

# print("Pass 3")
for (i in c(1:(totsteps-1)))
{
    # generate proposal
    z = rnorm(n = 1, sd = 0.25)
    prop = mcmc[i] + z
	
    # calculate likelihood and posterior density
    thetavec = truethetavec
    thetavec[1] = prop
    propden = Rdtq2d(thetavec, gammavec, runner, c1=xchase, c2=ychase, h = myh, numsteps = myns, k = myk, yM = xylimit)
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

