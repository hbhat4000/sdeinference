rm(list = ls(all = TRUE))

# solve inverse problem for linear SDE
# dX_t = theta1 (theta2 - X_t) dt +  theta3 dW_t

thetavec = c(0.5, 1, 0.25)

h = 1e-6
littlet = 5e-1
bigt = 250

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = 5e5
stopifnot((nsteps == (nsaves*hilt)))

ntrials = 1
h12 = sqrt(h)
xtraj = matrix(0, nrow = ntrials, ncol = (nsaves + 1))

# initial condition centered at origin so both -ve and +ve initial conditions generated
xtraj[,1] = rnorm(n = ntrials) 

for (i in c(1:nsaves))
{
  # print loop counter
  print(i)
  flush.console()
  
  x = xtraj[,i]
  
  for (j in c(1:hilt))
    # x_{t+1} = x_t + (\theta_1)*(\theta_2 - x_t)*h + (h12)*(\theta_3)*(z)
    x = x + (thetavec[1]*(thetavec[2] - x))*(h) + (h12)*(thetavec[3])*(rnorm(n = ntrials))
  
  xtraj[,(i+1)] = x
}

# y_t = x_t + \epsilon_t
xtraj = xtraj + matrix(rnorm(n = ((nsaves+1)*ntrials), mean = 0, sd = 0.1), nrow = ntrials, ncol = (nsaves+1))

tvec = seq(from = 0, to = bigt, by = littlet)
xtraj = rbind(tvec, xtraj)
save(xtraj, file = 'fakedata.RData')