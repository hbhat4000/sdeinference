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
theta = init_theta

# define log prior which computes the sum of log normals (prior)
logprior <- function(z)
{
    return(sum(dnorm(x = z, mean = hmc_mu, sd = hmc_mass, log = TRUE)))
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

    objective = -(sum(log(mylik)) + prior)
    gradient = numeric(length = numparam)

    # the other slices of the cube have the gradients with respect to each parameter
    for (i in c(1:numparam))
        gradient[i] = sum(probcube[,,(i+1)] / mylik)

    # add gradient of log prior to get gradient of log posterior
    gradient = gradient - (z - hmc_mu)/(hmc_mass^2)
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
    phi = rnorm(n = numparam, mean = hmc_mu, sd = hmc_mass)

    oldprior = logprior(phi)
    oldgradpost = logobjgrad(theta, phi, oldprior)

    # half step for momentum phi
    phi = phi + (oldgradpost$gradient) * (hmc_epsilon / 2)

    # L-1 full steps alternatively
    for(j in c(1:hmc_steps))
    {
        # full position step
        theta = theta + phi * (hmc_epsilon / hmc_mass)
        theta[1] = init_theta[1]

        # full momentum step except at end
        if(j != hmc_steps)
        {                      
            halfprior = logprior(phi)  
            halfgradpost = logobjgrad(theta, phi, halfprior)
            phi = phi + (halfgradpost$gradient) * (hmc_epsilon)
        }
    }

    propprior = logprior(phi)
    propgradpost = logobjgrad(theta, phi, propprior)
    intermediateval = propgradpost$objective + propprior - (oldgradpost$objective + oldprior)

    if (propgradpost$objective == -Inf)
        rho = 0
    else if (intermediateval > 0)
        rho = 1
    else
        rho = exp(intermediateval)

    # accept/reject step
    u = runif(n = 1)
    if (rho > u)
    {
        print(paste("Accepted step", i, ": ", paste("theta[", c(1:numparam), "]=", format(theta, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i] = 1
    }
    else
    {
        print(paste("Rejected step", i, ": ", paste("theta[", c(1:numparam), "]=", format(theta, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        artrack[i] = 0
    }
    thetamat[i,] = theta
    rmserror[i] = rmse(actualtheta, theta)
}

looptime = proc.time() - ptm

acceptance = 1 - mean(duplicated(thetamat[-(1:burnin),]))

########## First plotting technique ##########
# par(mfrow = c(2,3))
par(mfrow = c(2,2))

# hist(thetamat[-(1:burnin),1], nclass = 30, main = "Posterior of theta1", xlab = "True value = red line", xlim = c(actualtheta[1] - 0.05, actualtheta[1] + 0.05))
# abline(v = mean(thetamat[-(1:burnin),1]), col = "green")
# abline(v = actualtheta[1], col = "red")
# lines(density(thetamat[-(1:burnin),2]), col = "blue")

hist(thetamat[-(1:burnin),2], nclass = 30, main = "Posterior of theta2", xlab = "True value = red line", xlim = c(actualtheta[2] - 0.05, actualtheta[2] + 0.05))
abline(v = mean(thetamat[-(1:burnin),2]), col = "green")
abline(v = actualtheta[2], col = "red")
lines(density(thetamat[-(1:burnin),2]), col = "blue")

thetamat3new = (thetamat[-(1:burnin), 3])^2
newtheta3 = (actualtheta[3])^2
hist(thetamat3new, nclass = 30, main = "Posterior of theta3", xlab = "True value = red line", xlim = c(newtheta3 - 0.05, newtheta3 + 0.05))
abline(v = mean(thetamat3new), col = "green")
abline(v = newtheta3, col = "red")
lines(density(thetamat3new), col = "blue")

# hist(thetamat[-(1:burnin),3], nclass = 30, main = "Posterior of theta3", xlab = "True value = red line", xlim = c(actualtheta[3] - 0.05, actualtheta[3] + 0.05))
# abline(v = mean(thetamat[-(1:burnin),3]), col = "green")
# abline(v = actualtheta[3], col = "red")
# lines(density(thetamat[-(1:burnin),3]), col = "blue")

# plot(thetamat[-(1:burnin),1], type = "l", xlab = "True value = red line", ylab = "theta2 values", main = "HMC values of theta1", ylim = c(actualtheta[1] - 0.05, actualtheta[1] + 0.05))
# abline(h = actualtheta[1], col = "red")

plot(thetamat[-(1:burnin),2], type = "l", xlab = "True value = red line", ylab = "theta2 values", main = "HMC values of theta2", ylim = c(actualtheta[2] - 0.05, actualtheta[2] + 0.05))
abline(h = actualtheta[2], col = "red")

plot(thetamat3new, type = "l", xlab = "True value = red line", ylab = "theta3 values", main = "HMC values of theta3", ylim = c(newtheta3 - 0.05, newtheta3 + 0.05))
abline(h = newtheta3, col = "red")

# plot(thetamat[-(1:burnin),3], type = "l", xlab = "True value = red line", ylab = "theta3 values", main = "HMC values of theta3", ylim = c(actualtheta[3] - 0.05, actualtheta[3] + 0.05))
# abline(h = actualtheta[3], col = "red")

# for comparison:
# summary(lm(y~x))

########## Second plotting technique ##########
# Using MCMC package for plotting
library('coda')
chain = mcmc(thetamat)
summary(chain)
plot(chain)

myout = list("gridh" = gridh, "gridk" = gridk, "gridM" = gridM, "actualtheta" = actualtheta, "burnin" = burnin, "totsteps" = totsteps, "hmc_epsilon" = hmc_epsilon, "hmc_mu" = hmc_mu, "hmc_mass" = hmc_mass, "hmc_steps" = hmc_steps, "thetamat" = thetamat, "artrack" = artrack, "rmserror" = rmserror, "looptime" = looptime)
fname = paste('HMCsamples_', hmcdatanum, '.RData', sep = '')
save(myout, file = fname)
