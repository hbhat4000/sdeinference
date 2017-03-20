rm(list=ls(all = TRUE))

init = 2
final = 3

theta = c(1, 0.5, 1)

h = 0.01
k = 0.01
M = 1
deltat = 1
numsteps = ceiling(deltat/h)
grid = c(-M:M)*k

f <- function(theta, x) {
  return((theta[1])*(theta[2] - x))
}

g <- function(theta, x) {
  return(theta[3])
}

transition <- function(theta, x, y, t) {
  p = sqrt(theta[1] / (pi*(theta[3]^2)*(1 - exp(-2*theta[1]*t)))) * exp((-theta[1]*(x - theta[2] - (y - theta[2])*(exp(-theta[1]*t)))^2) / (theta[3]^2*(1 - exp(-2*theta[1]*t))))

  # specifically for \mu = 0
  # p = sqrt(theta[1] / (pi*(theta[3]^2*(1 - exp(-2*theta[1]*t))))) * exp((-theta[1]*(x - y*exp(-theta[1]*t))^2) / (theta[3]^2*(1 - exp(-2*theta[1]*t))))
  # print(p)
  return(p)
}

tdens = numeric(numsteps)

for(i in c(1:numsteps)) {
  tdens[i] = transition(theta, grid, init, i*h)
}

plot(tdens)
