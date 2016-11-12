# clear all memory
rm(list=ls(all=TRUE))

# load 2d dtq package
library('Rdtq2d')

# load fake data in X which has x1, x2, t
# load('fakedata_nonlinear.RData')

load('newfakedata_fullres.RData')
X = X[seq(from=1,to=2001,by=100),]

# keep every 100th row
# X = X[seq(from=1,to=nrow(X),by=100),]

# load true thetavec
source('truethetavec.R')

set.seed(1)

ptm = proc.time()

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
    loglik = sum(log(den[den >= 2.2e-16]))
    return(loglik + sum(log(prior)))
}

thetavec = truethetavec
thetavec[1] = x[1,1]
thetavec[2] = x[1,2]
oldden = Rdtq2d(thetavec,C1,C2,h=myh,numsteps=myns,k=myk,yM=xylimit)
checkk = PDFcheck(thetavec,h=myh,k=myk,yM=xylimit)
print(checkk)

oldpost = myposterior(den=oldden, prior=myprior(x[1,]))
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
    proppost = myposterior(den=propden, prior=myprior(prop)) 
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

print(proc.time() - ptm)


# percentage difference between empirical and true means
diffmeans = (mean(x[,1]^2) - 2*pi)/(2*pi)
print(diff)

# percentage difference between empirical and true modes
myden = density(x[,1]^2)
diffmodes = (myden$x[which.max(myden$y)] - 2*pi)/(2*pi)
print(diffmodes)

# save everything
save.image(file = 'samples1_by20_h_005.RData')
