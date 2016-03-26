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

objgradfun <- function(c0, z, mu, mass_sd)
{
    probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)

    mylik = probmat$lik
    mylik[mylik < 0] = 0

    objective = - sum(log(mylik))

    nc0 = length(c0)
    gradient = numeric(length = nc0)

    for (i in c(1:nc0))
        gradient[i] = - sum(probmat$grad[[i]] / probmat$lik)

    gradient = gradient - (z - mu)/(mass_sd)^2
    return(list("objective" = objective, "gradient" = gradient))
}

# create grid, compute densities on that grid
myh = 0.1
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.5))

# initial condition fakedata = c(1,4,0.5)
theta = c(1, 2, 1)
numparam = length(theta)

totsteps = 10
artrack = {}

for (i in c(1:(totsteps-1)))
{
    # generate proposal for the momentum term, phi is numparam dimensional 
    mass_sd = 0.25
    mu = 0
    phi = rnorm(n = numparam, mean = mu, sd = mass_sd)

    oldlogprior = logprior(phi, mu, mass_sd)
    oldgradpost = objgradfun(theta, phi, mu, mass_sd)

    phi_half = phi + (oldgradpost$gradient) * (myh / 2)
    theta_star = theta + phi_half * (myh / mass_sd)

    propgradpost = objgradfun(theta_star, phi_half, mu, mass_sd)
    phi_star = phi_half + (propgradpost$gradient) * (myh / 2)

    proplogprior = logprior(phi_star, mu, mass_sd)

	rho = (exp(propgradpost$objective) * exp(proplogprior)) / (exp(oldgradpost$objective) * exp(oldlogprior))

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        theta = theta_star
        oldgradpost = propgradpost
        print(paste("Accepted the proposal",prop))
        artrack[i] = 1
    }
    else
    {
        print(paste("Rejected the proposal",prop))
        artrack[i] = 0
    }
}