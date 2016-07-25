# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data which has 2 lists, chaser and runner, each with t, x, y
load('fakedata4.RData')
tchase = chaser[[1]]
xchase = chaser[[2]]
ychase = chaser[[3]]
runner = matrix(unlist(runner),ncol=3)

# algorithm parameters
mydatapoints = length(tchase) - 1
myh = 0.4
myk = 0.8*(myh)^0.9
xylimit = 5    # court dimension is 94*50
# xylimit = max(abs(xchase), abs(ychase))

xchase = xchase[1:mydatapoints]
ychase = ychase[1:mydatapoints]

# time increment from data
timeinc = tchase[2] - tchase[1]
myns = floor(timeinc/myh)

M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)
mm = length(xvec)

myposterior <- function(likden, dat, prior)
{
    steps = ncol(likden)
    datamat = matrix(c(dat[[2]], dat[[3]]), ncol = 2)
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
    logpost = logpost
    return(logpost)
}

posterior = {}
gamma1 = {}
for(i in seq(from = 1, to = 200, by = 1))
{ 
  gamma1[i] = i/100
  gammavec = c(gamma1[i], 1)
  nuvec = c(0.5, 0.5)

  oldden = Rdtq2d(nuvec, gammavec, runner, c1=xchase, c2=ychase, h = myh, numsteps = myns, k = myk, yM = xylimit)
  posterior[i] = myposterior(likden = oldden, dat = chaser, prior = myprior(c(gammavec, nuvec)))
  print(posterior[i])
}

plot(gamma1, posterior)

