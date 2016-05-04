rm(list=ls(all=TRUE))

library('nloptr')
load('fakedata.RData')
fd = xtraj
source('adjoint.R')

# define log prior which computes the sum of log normals (prior)
logprior <- function(z)
{
    return(sum(dnorm(x = z, mean = hmc_mu, sd = hmc_sd, log = TRUE)))
}

logobjgrad <- function(theta, z, prior)
{		
	out = gradobjfun(c0 = theta, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd,  diffcoeff = 1)

	# objective function returned with a negative for optimization
	myobj = out$obj
	mygrad = out$grad

	objective = sum(log(myobj)) + prior
    gradient = numeric(length = numparam)

    # add gradient of log prior to get gradient of log posterior
    gradient = mygrad - (z - hmc_mu)/(hmc_sd^2)

	return(list("objective" = objective, "gradient" = gradient))
}

myh = 0.1
myk = (myh)^(0.6)
mybigm = ceiling(pi/(myk^1.5))
theta = rnorm(12)

# samples
burnin = 200
samplesteps = 500
totsteps = burnin + samplesteps
numparam = length(theta)

# res <- nloptr(x0 = theta, eval_f = objgradfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level"=3, "check_derivatives" = FALSE, "xtol_abs"=1e-3))


thetamat = matrix(nrow = totsteps, ncol = numparam)
thetagenerated = matrix(nrow = totsteps, ncol = numparam)
artrack = numeric(length = totsteps)

##### HMC specific parameters #####
hh = 0.01
hmc_sd = 10
hmc_mu = 0

ptm = proc.time()
for (i in c(1:totsteps))
{
    # generate proposal for the momentum term, phi is numparam dimensional 
    phi = rnorm(n = numparam, mean = hmc_mu, sd = hmc_sd)

    oldlogprior = logprior(phi)
    oldgradpost = logobjgrad(theta, phi, oldlogprior)

    phi_half = phi + (oldgradpost$gradient) * (hh / 2)
    theta_star = theta + phi_half * (hh / hmc_sd)

    halflogprior = logprior(phi_half)
    propgradpost = logobjgrad(theta_star, phi_half, halflogprior)

    if (propgradpost$objective == -Inf)
        rho = 0
    else
    {
        phi_star = phi_half + (propgradpost$gradient) * (hh / 2)
        proplogprior = logprior(phi_star)
        rho = exp(propgradpost$objective + proplogprior - (oldgradpost$objective + oldlogprior))
    }

    thetagenerated[i,] = theta_star
    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        theta = theta_star
        oldgradpost = propgradpost
        print(paste("Accepted step", i, ": (", paste(format(theta_star, digits = 3, scientific = TRUE), collapse = ', ', sep = '')), ")")
        artrack[i] = 1
    }
    else
    {
        print(paste("Rejected step", i, ": ", paste(format(theta_star, digits = 3, scientific = TRUE), collapse = ', ', sep = '')), ")")
        artrack[i] = 0
    }
    thetamat[i,] = theta
}

looptime = proc.time() - ptm

acceptance = 1 - mean(duplicated(thetamat[-(1:burnin),]))

myout = list(myh, myk, bigm, theta, burnin, totsteps, hh, hmc_mu, hmc_sd, thetamat, artrack, rmserror, looptime)
fname = paste('Adjointsamples_', hmcdatanum, '.RData', sep = '')
save(myout, file = fname)