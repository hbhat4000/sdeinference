rm(list = ls(all = TRUE))

library(MASS)

# load data
load('fakedata_ou_noise.RData')
fd = xtraj

# load necessary functions
source('dtq_with_grad.R')

objgradfun <- function(c0)
{
    probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)
    mylik = probmat$lik
    mylik[mylik < 0] = 0
    objective = -sum(log(mylik))
    nc0 = length(c0)
    gradient = numeric(length=nc0)
    for (i in c(1:nc0))
        gradient[i] = -sum(probmat$grad[[i]] / probmat$lik)

    return(list("objective"=objective,"gradient"=gradient))
}

# create grid, compute densities on that grid
myh = 0.1
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.5))

# initial condition fakedata = c(1,4,0.5)
theta = c(1, 2, 1)
numparam = length(theta)

totsteps = 1000
# OLD: Metropolis, NEW: HMC

for (i in c(1:(totsteps-1)))
{
    # generate proposal for the momentum term, phi is numparam dimensional 
    # mu is a vector of means and sigma is the covariance matrix
    phi = mvrnorm(n = numparam, mu = 0, sd = 0.25)
	
#     # calculate likelihood and posterior density
#     thetavec = truethetavec
#     thetavec[1] = prop
#     propden = Rdtq2d(thetavec,C1,C2,h=myh,numsteps=myns,k=myk,yM=xylimit)
#     proppost = myposterior(likden=propden, dat=X, prior=myprior(prop)) 
#     rho = exp(proppost-oldpost)
#     maxcolsumerr = max(abs(colSums(propden)*myk^2 - 1))

#     # accept/reject step
#     u = runif(n = 1)
#     if (rho > u)
#     {
#         x[i+1] = prop
#         oldden = propden
#         oldpost = proppost
#         print(paste("Accepted the proposal",prop))
#         artrack[i] = 1
#     }
#     else
#     {
#         x[i+1] = x[i]
#         print(paste("Rejected the proposal",prop))
#         artrack[i] = 0
#     }
#     print(c(i,rho,maxcolsumerr,x[i+1]))
#     flush.console()
}