rm(list = ls(all = TRUE))

# load data
load('doublewell.RData')

# load necessary functions
source('dtq.R')

myh = 0.05
myk = myh^(0.75)
mybigm = ceiling(pi/(myk^1.5))
theta = c(1, 2, 0.7)
probmat = cdt(theta, h = myh, k = myk, bigm = mybigm, littlet = 1, data = xtraj)
mylik = probmat$lik
mylik[mylik < 0] = 0
objective = -sum(log(mylik))

