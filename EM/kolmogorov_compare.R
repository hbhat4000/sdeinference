rm(list=ls(all = TRUE))

init = 2
final = 3

theta = c(1, 1, 1)

h = 0.1
k = 0.01
M = 10
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

# T-t instead of t for backward
transition_backward <- function(theta, x, y, t, T) {
  p = sqrt(theta[1] / (pi*(theta[3]^2)*(1 - exp(-2*theta[1]*(T - t))))) * exp((-theta[1]*(x - theta[2] - (y - theta[2])*(exp(-theta[1]*(T-t))))^2) / (theta[3]^2*(1 - exp(-2*theta[1]*(T - t)))))
  return (p)
  }

tdens = numeric(length(grid))
savepdf = numeric(numsteps)

tdens = transition(theta, grid, init, h)
# plot(tdens)
savepdf[1] = tdens[M + 1]

for(i in c(2:numsteps - 1)) {
  tdens = transition(theta, grid, grid, i*h)
  # plot(tdens)
  savepdf[i] = tdens[M+1]
}

tdens = transition(theta, grid, final, numsteps*h)
# plot(tdens)
savepdf[numsteps] = tdens[M + 1]

print(savepdf)
plot(savepdf)
