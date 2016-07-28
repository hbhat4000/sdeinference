# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data which has 2 lists, chaser and runner, each with t, x, y
chaserlist = list(NULL)
runnerlist = list(NULL)
chasernew = list(NULL)
for (iii in c(1:100))
{
  fname = paste('./data/fakedata_h_0.04_',iii,'.RData',sep='')
  load(fname)
  chaserlist[[iii]] = matrix(unlist(chaser), ncol = 3)
  runnerlist[[iii]] = matrix(unlist(runner), ncol = 3)
  mydatapoints = nrow(chaserlist[[iii]]) - 1
  chasernew[[iii]] = (chaserlist[[iii]])[1:mydatapoints,]
  if (class(chasernew[[iii]])=="numeric") chasernew[[iii]] = t(as.matrix(chasernew[[iii]]))
}
print(max(unlist(chaserlist)))
print(max(unlist(runnerlist)))

# algorithm parameters
myh = 0.2
myk = 0.1 # 0.2*(myh)^0.9
xylimit = 22  # court dimension is 94*50

# time increment from data
timeinc = (runnerlist[[1]])[2,1] - (runnerlist[[1]])[1,1]
myns = floor(timeinc/myh)
print(myns)
M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)
mm = length(xvec)

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

lik = numeric(length=200)
gammavec = c(0,0.5)

for(i in seq(from = 5, to = 50, by = 5))
{ 
<<<<<<< HEAD:pursuit2d/demo2_local.R
<<<<<<< HEAD
  gammavec[1] = i/100
=======
# gamma1[i] = i/100
  gammavec = rep(i/100,2)
>>>>>>> 33550ec4e10fb8f972a9b5ea79a10de0007919b1
=======
# gamma1[i] = i/100
  gammavec = rep(i/100,2)
>>>>>>> 33550ec4e10fb8f972a9b5ea79a10de0007919b1:pursuit2d/demo2.R
  nuvec = c(0.5,0.5)

  for (iii in c(1:50))
  {
    oldden = Rdtq2d(nuvec, gammavec, runnerlist[[iii]], chasernew[[iii]], myh, myns, myk, xylimit)
    lik[i] = lik[i] + mylik(likden = oldden, dat = chaserlist[[iii]])
  }
  print(c(i/100,lik[i]))
}

plot(gamma1, lik)

