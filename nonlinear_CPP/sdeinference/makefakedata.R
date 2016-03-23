rm(list=ls(all=TRUE))

f <- function(X,thetavec) 
{
  return(0.5*(thetavec[1] - X))
}

g <- function(X,thetavec)
{
  return(thetavec[2])
}

# Compute samples of 2-D SDE
# dX_t =  0.5*(theta1 - X_t) dt + theta2 dW_t
# at different time, with final time T = 1

source('truethetavec.R')
dt = 0.001
T = 1
nsteps = ceiling(T/dt) + 1

# columns are x, t
X = matrix(0, nrow = nsteps, ncol = 2)
X[1,] = c(0.1,0)

for (i in c(2:(nsteps)))
{
  z = rnorm(n = 1, mean = 0, sd = sqrt(dt))
  X[i,1] = X[i-1,1] + f(X[i-1,1],truethetavec)*dt + g(X[i-1,1],truethetavec)*z
  X[i,2] = X[i-1,2] + dt
}

X = X[seq(from = 1, to = nsteps, by = 25),]
save(X, file = 'fakedata.RData')