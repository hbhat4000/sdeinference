rm(list = ls(all = TRUE))

# load data
load('fakedata.RData')
fd = xtraj[1:4,]

# load necessary functions
source('adjoint.R')

myh = 0.01
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.5))

# initial condition fakedata = c(1,4,0.5)
theta = c(1, 0.5, 0.5)

# run the adjoint code
test = gradobjfun(theta, myh, myk, mybigm, 1, fd)
print(test)

