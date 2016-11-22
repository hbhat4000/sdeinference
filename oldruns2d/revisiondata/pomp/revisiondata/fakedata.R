rm(list=ls(all=TRUE))

f1 <- function(X,thetavec) 
{
  return(-X[2]*thetavec[1])
}
f2 <- function(X,thetavec) 
{
  return(X[1]*thetavec[2])
}
g1 <- function(X,thetavec)
{
  return((thetavec[3]^2)*thetavec[1])
}
g2 <- function(X,thetavec)
{
  return((thetavec[4]^2)*thetavec[2])
}

# Compute samples of 2-D SDE
# dX_t =  theta1 Y_t dt + theta3^2 dW_t
# dY_t =  theta2 X_t dt + theta4^2 dW_t
# at different time, with final time T = 4

source('truethetavec.R')
dt = 0.001
T = 2
nsteps = ceiling(T/dt) + 1

# columns are x1, x2, t
X = matrix(0, nrow = nsteps, ncol = 3)
X[1,] = c(.1,0,0)

for (i in c(2:(nsteps)))
{
  z1 = rnorm(n = 1, mean = 0, sd = sqrt(dt))
  z2 = rnorm(n = 1, mean = 0, sd = sqrt(dt))
  X[i,1] = X[i-1,1] + f1(X[i-1,],truethetavec)*dt + g1(X[i-1,],truethetavec)*z1
  X[i,2] = X[i-1,2] + f2(X[i-1,],truethetavec)*dt + g2(X[i-1,],truethetavec)*z2
  X[i,3] = X[i-1,3] + dt
}

save(X, file='newfakedata_fullres.RData')


