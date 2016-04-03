# Compute samples of 2-D SDE
# dX_t =  -theta1 Y_t dt + (theta3)^2 dW_t
# dY_t =  theta2 X_t dt + (theta4)^2 dW_t

source('truethetavec.R')

# f1 <- function(x1, x2, thetavec) 
# {
#   return(-thetavec[1]*x2)
# }

# f2 <- function(x1, x2, thetavec) 
# {
#   return(thetavec[2]*x1)
# }

# g1 <- function(x1, x2, thetavec)
# {
#   return(thetavec[3]^2)
# }

# g2 <- function(x1, x2, thetavec)
# {
#   return(thetavec[4]^2)
# }

f1 <- function(x, thetavec) 
{
  return(-thetavec[1]*x[2])
}

f2 <- function(x, thetavec) 
{
  return(thetavec[2]*x[1])
}

g1 <- function(x, thetavec)
{
  return(thetavec[3]^2)
}

g2 <- function(x, thetavec)
{
  return(thetavec[4]^2)
}

driftfun <- function(thetavec, y) 
{
  # computes drift as a matrix
  f = ((thetavec[1])*(y)*(thetavec[2]) - y^2)
  return(f)
}

difffun <- function(thetavec, y)
{
  # replicates diffusion as a long vector
  g = rep(thetavec[3], length(y))
  return(g)
}

difffun1 <- function(thetavec, y)
{
  # replicates diffusion term as a matrix
  dd1 = dim(y)[1]
  dd2 = dim(y)[2]
  vec = rep(thetavec[3], dd1)
  return(replicate(dd2, vec)) 
}
