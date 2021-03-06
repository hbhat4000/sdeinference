rm(list = ls(all = TRUE))

# load data
load('fakedata_ou_noise.RData')
fd = xtraj

# load necessary functions
source('dtq_with_grad.R')

# define log prior
logprior <- function(z, mu, mass_sd)
{
    return(dnorm(x = z, mean = mu, sd = mass_sd, log = TRUE))
}

objgradfun <- function(c0, z, mu, mass_sd, prior)
{
    probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)

    mylik = probmat$lik
    mylik[mylik < 0] = 0

    # log posterior = log likelihood + log prior
    objective = sum(log(mylik)) + prior

    # form gradient of log likelihood
    nc0 = length(c0)
    gradient = numeric(length = nc0)
    for (i in c(1:nc0))
        gradient[i] = sum(probmat$grad[[i]] / probmat$lik)

    # add gradient of log prior to get gradient of log posterior
    gradient = gradient - (z - mu)/(mass_sd^2)
    return(list("objective" = objective, "gradient" = gradient))
}

# create grid, compute densities on that grid
myh = 0.01
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.5))

# initial condition fakedata = c(1,4,0.5)
theta = c(1, 2, 1)
numparam = length(theta)

hh = 0.005
totsteps = 1000
thetamat = matrix(nrow=totsteps, ncol=numparam)
artrack = numeric(length=totsteps)

for (i in c(1:totsteps))
{
    # generate proposal for the momentum term, phi is numparam dimensional 
    mass_sd = 2
    mu = 0
    phi = rnorm(n = numparam, mean = mu, sd = mass_sd)

    oldlogprior = sum(logprior(phi, mu, mass_sd))
    oldgradpost = objgradfun(theta, phi, mu, mass_sd, oldlogprior)

    phi_half = phi + (oldgradpost$gradient) * (hh / 2)
    theta_star = theta + phi_half * (hh / mass_sd)

    halflogprior = sum(logprior(phi_half, mu, mass_sd))
    propgradpost = objgradfun(theta_star, phi_half, mu, mass_sd, halflogprior)

    if (propgradpost$objective == -Inf)
        rho = 0
    else
    {
        phi_star = phi_half + (propgradpost$gradient) * (hh / 2)
        proplogprior = sum(logprior(phi_star, mu, mass_sd))
        rho = exp(propgradpost$objective + proplogprior - (oldgradpost$objective + oldlogprior))
    }

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        theta = theta_star
        oldgradpost = propgradpost
        print(paste("Accepted",paste("theta[",c(1:3),"]=",format(theta_star,digits=3,scientific=TRUE),collapse=', ',sep='')))
        artrack[i] = 1
    }
    else
    {
        print(paste("Rejected",paste("theta[",c(1:3),"]=",format(theta_star,digits=3,scientific=TRUE),collapse=', ',sep='')))
        artrack[i] = 0
    }
    thetamat[i,] = theta
}
myout = list(theta=thetamat,ar=artrack)
fname = paste('posteriorsamples_',myh,'.RData',sep='')
save(myout,file=fname)

