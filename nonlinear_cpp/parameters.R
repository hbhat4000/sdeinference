# index of the data set being used (to maintain consistency)
datanum = 1

# parameters for fakedata
fakedatah = 1e-4
littlet = 1
bigt = 50
ntrials = 100

# true parameter for RMSE calculation
actualtheta = c(1, 4, 0.5)

# parameters for inference
# create grid, compute densities on that grid
gridh = 0.05
gridk = (gridh)^0.85
gridM = ceiling(pi/(gridk)^1.5)

# initial condition
theta = c(2, 2, 2)
numparam = length(theta)

# samples
burnin = 100
samplesteps = 1000
totsteps = burnin + samplesteps

# MCMC specific parameters
prop_mu = 0
prop_sd = 0.1
prior_mu = 0
prior_sd = 100

# HMC specific parameters
hh = 0.01
hmc_sd = 2
hmc_mu = 0