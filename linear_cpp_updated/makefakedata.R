rm(list = ls(all = TRUE))

# load all required parameters (initial theta, grid size, ...)
source('parameters.R')

# solve inverse problem for nonlinear SDE
# dX_t = theta1 X_t (theta2 - (X_t)^2) dt +  theta3 dW_t

# theta1, theta2 > 0
# stable equilibrium at +sqrt(theta2) or -sqrt(theta2) depending on the IC
thetavec = actualtheta
h = fakedatah

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

h12 = sqrt(h)
xtraj = matrix(0, nrow = ntrials, ncol = (nsaves + 1))

# initial condition centered at origin so both -ve and +ve initial conditions generated
xtraj[,1] = runif(n = ntrials, min = 0, max = 2) 

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x = xtraj[,i]

    for (j in c(1:hilt))
        x = x + (thetavec[1]*(thetavec[2] - x))*(h) + (h12)*(thetavec[3]*thetavec[3])*(rnorm(n = ntrials))

    xtraj[,(i+1)] = x
}

# tvec = seq(from = 0, to = bigt, by = littlet)
# xtraj = rbind(tvec, xtraj)

fname = paste('fakedata_', fakedatanum, '.RData', sep = '')
save(xtraj, file = fname)
