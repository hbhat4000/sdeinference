fakedatanum = 4
mcmcdatanum = 4
hmcdatanum = 4

# parameters for fakedata
fakedatah = 1e-4
littlet = 1
bigt = 25
ntrials = 50

# true parameter for RMSE calculation
actualtheta = c(0.8, 0.9, 0.5)

# parameters for inference
# create grid, compute densities on that grid=
gridh = 0.05
gridk = (gridh)^0.75
gridM = ceiling(pi/(gridk)^1.5)

# initial condition
theta = c(1, 2, 1)
numparam = length(theta)

# samples
burnin = 200
samplesteps = 500
totsteps = burnin + samplesteps

##### MCMC specific parameters #####
prop_mu = 0
# prop_sd = 0.05
prop_sd = c(0.01, 0.05, 0.05)
prior_mu = 0
prior_sd = 100

##### HMC specific parameters #####
hh = 0.01
hmc_sd = 10
hmc_mu = 0

######## Comments ##########
# with prop_sd = 0.001, a/r ratio = 147.1% (655/445)
# with prop_sd = 0.05, a/r ratio =  9.34%
# with prop_sd = 0.01, a/r ratio = 4.56%


####### Locally installing directories #########
# [shagun@ucmerced-169-236-72-226 linear_cpp]$ export R_LIBS="/home/shagun/R"
# [shagun@ucmerced-169-236-72-226 linear_cpp]$ R CMD INSTALL -l /home/shagun/R Rgdtq
