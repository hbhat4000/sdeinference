fakedatanum = 5
mcmcdatanum = 5
hmcdatanum = 5

# parameters for fakedata
fakedatah = 1e-6
littlet = 1
bigt = 25
ntrials = 100

# true parameter for RMSE calculation
actualtheta = c(2, 3, 0.5)

# parameters for inference
# create grid, compute densities on that grid=
gridh = 0.05
gridk = (gridh)^0.75
gridM = ceiling(pi/(gridk)^1.5)

# initial condition
theta = c(1, 1, 1)
numparam = length(theta)

# samples
burnin = 1000
samplesteps = 2000
totsteps = burnin + samplesteps

##### MCMC specific parameters #####
prop_mu = 0
# prop_sd = 0.05
prop_sd = c(0.01, 0.01, 0.01)
prior_mu = 0
prior_sd = 100

##### HMC specific parameters #####
hh = 0.01
hmc_sd = 2
hmc_mu = 0

######## Comments ##########
# with prop_sd = 0.001, a/r ratio = 147.1% (655/445)
# with prop_sd = 0.05, a/r ratio =  9.34%
# with prop_sd = 0.01, a/r ratio = 4.56%
