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
  fname = paste('fakedata100/fakedata_h_0.4_',iii,'.RData',sep='')
  load(fname)
  chaserlist[[iii]] = matrix(unlist(chaser), ncol = 3)
  runnerlist[[iii]] = matrix(unlist(runner), ncol = 3)
  mydatapoints = nrow(chaserlist[[iii]]) - 1
  chasernew[[iii]] = (chaserlist[[iii]])[1:mydatapoints,]
  if (class(chasernew[[iii]])=="numeric") chasernew[[iii]] = t(as.matrix(chasernew[[iii]]))
}

# algorithm parameters
myh = 0.4
myk = 0.05 # 0.05*(myh)^(1.2)
xylimit = ceiling(max(unlist(chaserlist), unlist(runnerlist)))

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


lik = numeric(length = 200)
gammavec = c(0,0.5)

for(i in seq(from = 10, to = 85, by = 1))
{
  gammavec = rep(i/100,2)
  print(gammavec)
  nuvec = c(0.5,0.5)

  # gammavec = c(0.2, 0.5)
  # nuvec = rep(i/100, 2)
  # # print(nuvec)

  for (iii in c(1:100))
  {
    oldden = Rdtq2d(nuvec, gammavec, runnerlist[[iii]], chasernew[[iii]], myh, myns, myk, xylimit)
    lik[i] = lik[i] + mylik(likden = oldden, dat = chaserlist[[iii]])
  }
  print(c(i/100,lik[i]))
}

print(which.max(lik)/100)

plot(seq(from = 10, to = 85, by = 1)/100, lik)
