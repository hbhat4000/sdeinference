# rm(list=ls(all = TRUE))
# 
# init = 1
# final = 3
# 
# theta = c(1, 0.5, 2)
# 
# h = 0.1
# k = 0.01
# M = 800
# deltat = 1
# numsteps = ceiling(deltat/h)
grid = c((-M):M)*k

f <- function(theta, x) {
  return((theta[1])*(theta[2] - x))
}

g <- function(theta, x) {
  return(theta[3])
}

# keep in mind that y is the initial point
# and x is the final point
transition_forward <- function(theta, x, y, t) {
  p = sqrt(theta[1] / (pi*(theta[3]^2)*(1 - exp(-2*theta[1]*t)))) * exp((-theta[1]*(x - theta[2] - (y - theta[2])*(exp(-theta[1]*t)))^2) / (theta[3]^2*(1 - exp(-2*theta[1]*t))))

  # specifically for \mu = 0
  # p = sqrt(theta[1] / (pi*(theta[3]^2*(1 - exp(-2*theta[1]*t))))) * exp((-theta[1]*(x - y*exp(-theta[1]*t))^2) / (theta[3]^2*(1 - exp(-2*theta[1]*t))))
  # print(p)
  return(p)
}

transition_backward <- function(theta, x, y, t) {
  p = sqrt(theta[1] / (pi*(theta[3]^2)*(1 - exp(-2*theta[1]*t)))) * exp((-theta[1]*(x - theta[2] - (y - theta[2])*(exp(-theta[1]*t)))^2) / (theta[3]^2*(1 - exp(-2*theta[1]*t))))
  
  # specifically for \mu = 0
  # p = sqrt(theta[1] / (pi*(theta[3]^2*(1 - exp(-2*theta[1]*t))))) * exp((-theta[1]*(x - y*exp(-theta[1]*t))^2) / (theta[3]^2*(1 - exp(-2*theta[1]*t))))
  # print(p)
  print(p)
  return(p)
}

# tdens = numeric(2*M+1)
# savepdf = numeric(numsteps)
# 
# for(i in c(1:numsteps)) {
#   timestep = i*h
#   tdens = transition_forward(theta, grid, init, timestep)
#   savepdf[i] = tdens[M+1]
# }
# 
# plot(savepdf)
