# clear all memory
rm(list=ls(all=TRUE))

# load fake data in X which has x1, x2, t
load('./vanderpol/fakedata_vanderpol_fullres_new.RData')
mydata = X[seq(from = 1, to = nrow(X), by = 100),]

set.seed(1)

ptm = proc.time()

# algorithm parameters
# time increment from data and time step
timeinc = mydata[2,3] - mydata[1,3]
myh = timeinc/4

numparticles = 1000
source('pompinit.R')

burnin = 100
numsteps = 25000
totsteps = numsteps + burnin

mcmcsamples = matrix(0, nrow = totsteps, ncol = 3)
mcmcsamples[1,] = c(0.1,0.1,0.1)

# define log prior
myprior <- function(z)
{
    return(dnorm(x = z, mean = 0, sd = 100, log = TRUE))
}

# function that computes log posterior
myposterior <- function(loglik, prior)
{
    return(loglik + sum(prior))
}

oldpf <- pfilter(mymod, Np = numparticles, params = c(theta1 = mcmcsamples[1,1], theta2 = mcmcsamples[1,2], theta3 = mcmcsamples[1,3], X1.0 = mymod.dat[1,]$Y1, X2.0 = mymod.dat[1,]$Y2))
oldden <- logLik(oldpf)
oldpost = myposterior(oldden, prior=myprior(mcmcsamples[1,]))
artrack = numeric(length=(totsteps-1))

adaptivesd = 0.25

for (i in c(1:(totsteps-1)))
{
    z = rnorm(n=3, sd=adaptivesd)
    prop = mcmcsamples[i,] + z

    # calculate likelihood and posterior density
    proppf <- pfilter(mymod, Np = numparticles, params = c(theta1 = prop[1], theta2 = prop[2], theta3 = prop[3], X1.0 = mymod.dat[1,]$Y1, X2.0 = mymod.dat[1,]$Y2))
    propden <- logLik(proppf)
    proppost = myposterior(propden, prior = myprior(prop)) 
    # print(c(oldden,propden))
    # print(c(oldpost,proppost))
    rho = exp(proppost-oldpost)

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        mcmcsamples[i+1,] = prop
        oldden = propden
        oldpost = proppost
        # print(paste("Accepted step", i, ": ", paste("theta[", c(1:3), "]=", format(prop, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i+1] = 1
    }
    else
    {
        mcmcsamples[(i+1),] = mcmcsamples[i,]
        # print(paste("Rejected step", i, ": ", paste("theta[", c(1:3), "]=", format(prop, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i+1] = 0
    }
    if (i %% 100 == 0)
    {
      print(c(i,mean(artrack[1:i])))
    }
    flush.console()
}

finaltime = proc.time() - ptm
print(finaltime)

artrack = artrack[(burnin+1):totsteps]
arratio = sum(artrack)/length(artrack)

# throw away the burnin steps
mcmcsamples = mcmcsamples[(burnin+1):totsteps,]

# save everything
save.image(file = 'samples_pomp_noadapt_vanderpol_by4.RData')

