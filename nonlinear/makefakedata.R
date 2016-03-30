rm(list = ls(all = TRUE))

# solve inverse problem for nonlinear SDE
# dX_t = theta1 X_t (theta2 - (X_t)^2) dt +  theta3 dW_t

# theta1, theta2 > 0
# stable equilibrium at +sqrt(theta2) or -sqrt(theta2) depending on the IC
thetavec = c(1, 1, 1)

h = 0.0001
littlet = 1
bigt = 25

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

ntrials = 300
h12 = sqrt(h)
xtraj = matrix(0, nrow = ntrials, ncol = (nsaves + 1))

# initial condition centered at origin so both -ve and +ve initial conditions generated
xtraj[,1] = rnorm(n = ntrials, mean = 0, sd = 1) 

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x = xtraj[,i]

    for (j in c(1:hilt))
        x = x + (thetavec[1])*(x)*(thetavec[2] - x^2)*(h) + (h12)*(thetavec[3])*(rnorm(n = ntrials))      
    xtraj[,(i+1)] = x
}

tvec = seq(from = 0, to = bigt, by = littlet)
xtraj = rbind(tvec, xtraj)
save(xtraj, file = 'fakedata5.RData')

# Initial condition picked is printed out 
# print(xtraj[2,1])
