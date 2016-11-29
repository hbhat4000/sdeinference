# clear all memory
rm(list=ls(all=TRUE))

# load 2d dtq package
library('Rdtq2d')

# load fake data in X which has x1, x2, t
load('./vanderpol/fakedata_vanderpol_fullres_new.RData')
mydata = X[seq(from = 1, to = nrow(X), by = 100),]

# load true thetavec
source('./vanderpol/truethetavec_vanderpol.R')

set.seed(1)

ptm = proc.time()

# algorithm parameters
# time increment from data and time step
timeinc = mydata[2,3] - mydata[1,3]
myh = timeinc/16
myns = floor(timeinc/myh)
# print(myh)
# myk = 0.8*myh^0.75
myk = myh^(0.55)
xylimit = 1.5*max(abs(mydata[,1:2]))

xcord = mydata[,1]
ycord = mydata[,2]

burnin = 1000
numsteps = 20000
totsteps = numsteps + burnin

mcmcsamples = matrix(0, nrow = totsteps, ncol = 3)
mcmcsamples[1,] = c(0.1,0.1,0.1)
thetavec = truethetavec
thetavec[1:3] = mcmcsamples[1,]

M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)
mm = length(xvec)

# define log prior
myprior <- function(z)
{
    return(dnorm(x = z, mean = 0, sd = 100, log = TRUE))
}

# function that computes log posterior
myposterior <- function(den, prior)
{
	den[den <= 2.2e-16] = 2.2e-16
    loglik = sum(log(den))
    return(loglik + sum(prior))
}

oldden = Rdtq2d(thetavec, xcord, ycord, h = myh, numsteps = myns, k = myk, yM = xylimit)
# print(oldden)
# checkk = PDFcheck(thetavec,h=myh,k=myk,yM=xylimit)
# print("Did PDFcheck")
# print(checkk)

oldpost = myposterior(den=oldden, prior=myprior(mcmcsamples[1,]))
artrack = numeric(length=(totsteps-1))

for (i in c(1:(totsteps-1)))
{
    # generate proposal
    z = rnorm(n = 3, sd = 0.025)
    prop = mcmcsamples[i,] + z
	
    # calculate likelihood and posterior density
    # thetavec = truethetavec
    thetavec[1:3] = prop
    propden = Rdtq2d(thetavec, xcord, ycord, h = myh, numsteps = myns, k = myk, yM = xylimit)
    proppost = myposterior(den = propden, prior = myprior(prop)) 
    rho = exp(proppost-oldpost)
    maxcolsumerr = max(abs(colSums(propden)*myk^2 - 1))

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        mcmcsamples[i+1,] = prop
        oldden = propden
        oldpost = proppost
        # print(paste("Accepted step", i, ": ", paste("theta[", c(1:2), "]=", format(prop, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i+1] = 1
    }
    else
    {
        mcmcsamples[i+1,] = mcmcsamples[i,]
        # print(paste("Rejected step", i, ": ", paste("theta[", c(1:2), "]=", format(prop, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i+1] = 0
    }
    # print(c(i,rho,maxcolsumerr,x[i+1,]))
    if ((i %% 100) == 0) print(c(i,mean(artrack[1:i])))
    flush.console()
}

finaltime = proc.time() - ptm
print(finaltime)

artrack = artrack[(burnin+1):totsteps]
arratio = sum(artrack)/length(artrack)
# throw away the burnin steps
mcmcsamples = mcmcsamples[(burnin+1):totsteps,]

# percentage difference between empirical and true means
diffmeans = (mean(mcmcsamples[,1]^2) - 2*pi)/(2*pi)
print(diffmeans)

# percentage difference between empirical and true modes
myden = density(mcmcsamples[,1]^2)
diffmodes = (myden$x[which.max(myden$y)] - 2*pi)/(2*pi)
print(diffmodes)

# save everything
save.image(file = 'samples_vanderpol_by16.RData')

