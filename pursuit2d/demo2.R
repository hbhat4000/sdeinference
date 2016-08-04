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
  fname = paste('fakedataswitch/fakedata_h_0.2_',iii,'.RData',sep='')
  load(fname)
  chaserlist[[iii]] = matrix(unlist(chaser), ncol = 3)
  runnerlist[[iii]] = matrix(unlist(runner), ncol = 3)
  mydatapoints = nrow(chaserlist[[iii]]) - 1
  chasernew[[iii]] = (chaserlist[[iii]])[1:mydatapoints,]
  if (class(chasernew[[iii]])=="numeric") chasernew[[iii]] = t(as.matrix(chasernew[[iii]]))
}

# algorithm parameters
myh = 0.2
myk = 0.1 # 0.05*(myh)^(1.2)
xylimit = 1.5*ceiling(max(unlist(chaserlist), unlist(runnerlist)))

# time increment from data
timeinc = (runnerlist[[1]])[2,1] - (runnerlist[[1]])[1,1]
myns = floor(timeinc/myh)
M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)
mm = length(xvec)

# print(xylimit)
# print(timeinc)
# print(myns)
# print(M)
# print(mm)

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
#       likdat[likdat <= 2.2e-16] = 2.2e-16
      logpost = logpost + sum(log(likdat[likdat >= 2.2e-16]))
    }
    return(logpost)
}


lik = numeric(length = 200)

for(i in seq(from = 90, to = 140, by = 5))
{
  # gammavec = rep(i/100,2)
  # gammavec = c(i/100, 1.2)
  # gammavec = c(0.25, i/100)
  # gammavec = c(0.25, 0.75)
  # nuvec = c(i/100, 0.5)
  # gammavec = c(0.5,1.2)
  gammavec = c(0.5,i/100)
  nuvec = c(0.1,0.2)
  # print(nuvec)
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

liknew = lik[-which(lik==0)]
print(which.max(liknew))

# plot(seq(from = 10, to = 55, by = 5)/100, liknew)
