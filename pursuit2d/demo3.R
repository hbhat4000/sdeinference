# clear all memory
rm(list = ls(all = TRUE))

# load 2d dtq package
library('Rdtq2d')

# load package for interpolation
library('fields')

# load fake data which has 2 lists, chaser and runner, each with t, x, y
chaserlist = list(NULL)
runnerlist = list(NULL)

xymax = 0
for (iii in c(1:1))
{
  fname = paste('fakedataswitch4/fakedata_dt_0.1.RData',sep='')
  load(fname)
  chaserlist[[iii]] = matrix(unlist(chaser), ncol = 3)
  runnerlist[[iii]] = matrix(unlist(runner), ncol = 3)
  xymax = max(xymax,max(chaserlist[[iii]][,-1]))
  xymax = max(xymax,max(runnerlist[[iii]][,-1]))
}

# time increment from data
timeinc = (runnerlist[[1]])[2,1] - (runnerlist[[1]])[1,1]

# DTQ algorithm parameters
myh = timeinc/4
myk = timeinc/4 # 0.05*(myh)^(1.2)
xylimit = 1.5*xymax

# number of DTQ steps
myns = floor(timeinc/myh)

# M = ceiling(xylimit/myk)
# xvec = myk*c(-M:M)
# mm = length(xvec)

# print(xylimit)
# print(timeinc)
# print(myns)
# print(M)
# print(mm)

lik = numeric(length = 200)
gammalen = 8/myh

for(i in seq(from = 80, to = 120, by = 5))
{
  gammavec = c(rep(0.4,gammalen/2),rep(i/100,gammalen/2))
  nuvec = c(0.15,0.1)
  for (iii in c(1:1))
  {
    den = Rdtq2d(nuvec, gammavec, runnerlist[[iii]], chaserlist[[iii]], myh, myns, myk, xylimit)
    loglik = sum(log(den[den >= 2.2e-16]))
    lik[i] = lik[i] + loglik
  }
  print(c(i/100,lik[i]))
}

liknew = lik[-which(lik==0)]
print(which.max(liknew))

# plot(seq(from = 10, to = 55, by = 5)/100, liknew)
