# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data which has 2 lists, chaser and runner, each with t, x, y
load('fakedata5.RData')
chaser = matrix(unlist(chaser), ncol = 3)
runner = matrix(unlist(runner), ncol = 3)

# algorithm parameters
mydatapoints = nrow(chaser) - 1
myh = 0.4
myk = 0.8*(myh)^0.9
xylimit = 5    # court dimension is 94*50

# time increment from data
timeinc = runner[2,1] - runner[1,1]
myns = floor(timeinc/myh)
print(myns)
M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)
mm = length(xvec)

chasernew = chaser[1:mydatapoints, ]

mylik <- function(likden, dat)
{
    steps = ncol(likden)
    datamat = dat[,-1]
    logpost = 0

    for(i in c(1:steps))
    {
      myloc = matrix(datamat[(i+1),c(1:2)], ncol = 2)
      # print(myloc)

      # interp.surface(obj, loc)
      # obj : a list with components x, y and z s.t. x and y are the X and Y grid values and z is a matrix
      # with the corresponding values of the surface
      # loc : a matrix of irregular locations to interpolate, first column of loc is the X coordinates
      # and second column is the Y's
      likdat = interp.surface(obj = list(x = xvec, y = xvec, z = matrix(likden[,i], nrow = mm)), loc = myloc)
      likdat[likdat <= 2.2e-16] = 2.2e-16
      logpost = logpost + sum(log(likdat))
    }
    return(logpost)
}

lik = {}
gammavec = rep(1, myns)
gamma1 = {}

for(i in seq(from = 1, to = 200, by = 1))
{ 
  gamma1[i] = i/100
  gammavec[1] = gamma1[i]
  nuvec = c(0.5,0.5)

  oldden = Rdtq2d(nuvec, gammavec, runner, chasernew, myh, myns, myk, xylimit)
  lik[i] = mylik(likden = oldden, dat = chaser)
  print(lik[i])
}

plot(gamma1, lik)

