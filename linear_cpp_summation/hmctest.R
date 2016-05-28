rm(list = ls(all = TRUE))

# load all required parameters (initial theta, grid size, ...)
source('parameters.R')

# load the dtq package for finding the computing likelihood function
library('Rgdtq')

# library for RMSE calculation
library('Metrics')

# load data which is saved as fd 
fname = paste('fakedata_', fakedatanum, '.RData', sep = '')
load(fname)
fd = xtraj

# define log prior which computes the sum of log normals (prior)
logprior <- function(z)
{
    return(sum(dnorm(x = z, mean = hmc_mu, sd = hmc_sd, log = TRUE)))
}

# logposterior = log likelihood + log prior 
# likelihood is computed using the DTQ method by calling the Rgdtq C++ function
logobjgrad <- function(theta, z, prior)
{
    numparam = length(theta)
    # probmat = array(0, dim = c(nrow(fd), ncol(fd) - 1, numparam + 1))
    probcube = Rgdtql(thetavec = theta, h = gridh, k = gridk, M = gridM, littlet = 1, init_data = fd)

    # the first slice of cube is the likelihood
    mylik = probcube[,,1]
    mylik[mylik < 0] = 0

    objective = sum(log(mylik)) + prior
    gradient = numeric(length = numparam)

    # the other slices of the cube have the gradients with respect to each parameter
    for (i in c(1:numparam))
        gradient[i] = sum(probcube[,,(i+1)] / mylik)

    # add gradient of log prior to get gradient of log posterior
    gradient = gradient - (z - hmc_mu)/(hmc_sd^2)
    return(list("objective" = objective, "gradient" = gradient))
}


thetamat = matrix(nrow = totsteps, ncol = numparam)
thetagenerated = matrix(nrow = totsteps, ncol = numparam)
artrack = numeric(length = totsteps)
rmserror = numeric(length = totsteps)

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
        print(paste("Accepted step", i, ": ", paste("theta[", c(1:numparam), "]=", format(theta_star, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i] = 1
    }
    else
    {
        print(paste("Rejected step", i, ": ", paste("theta[", c(1:numparam), "]=", format(theta_star, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i] = 0
    }
    thetamat[i,] = theta
    rmserror[i] = rmse(actualtheta, theta)
}

looptime = proc.time() - ptm

acceptance = 1 - mean(duplicated(thetamat[-(1:burnin),]))

par(mfrow = c(2,3))

hist(thetamat[-(1:burnin),1], nclass = 30, main = "Posterior of theta1", xlab = "True value = red line")
abline(v = mean(thetamat[-(1:burnin),1]))
abline(v = actualtheta[1], col = "red")

hist(thetamat[-(1:burnin),2], nclass = 30, main = "Posterior of theta2", xlab = "True value = red line")
abline(v = mean(thetamat[-(1:burnin),2]))
abline(v = actualtheta[2], col = "red")

hist(thetamat[-(1:burnin),3], nclass = 30, main = "Posterior of theta3", xlab = "True value = red line")
abline(v = mean(thetamat[-(1:burnin),3]))
abline(v = actualtheta[3], col = "red")

plot(thetamat[-(1:burnin),1], type = "l", xlab = "True value = red line", main = "HMC values of theta1")
abline(h = actualtheta[1], col = "red")

plot(thetamat[-(1:burnin),2], type = "l", xlab = "True value = red line", main = "HMC values of theta2")
abline(h = actualtheta[2], col = "red")

plot(thetamat[-(1:burnin),3], type = "l", xlab = "True value = red line", main = "HMC values of theta3")
abline(h = actualtheta[3], col = "red")

# for comparison:
# summary(lm(y~x))

myout = list(gridh, gridk, gridM, actualtheta, burnin, totsteps, hh, hmc_mu, hmc_sd, thetamat, artrack, rmserror, looptime)
fname = paste('HMCsamples_', hmcdatanum, '.RData', sep = '')
save(myout, file = fname)
