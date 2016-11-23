# clear all memory
rm(list=ls(all=TRUE))

# load 2d dtq package
library('Rdtq2d')

# load fake data in X which has x1, x2, t
load('./vanderpol/fakedata_vanderpol_fullres_new.RData')
mydata = X[seq(from = 1, to = nrow(X), by = 10),]

# load true thetavec
source('./vanderpol/truethetavec_vanderpol.R')

set.seed(1)

# algorithm parameters
# time increment from data and time step
timeinc = mydata[2,3] - mydata[1,3]
myh = timeinc/1
myns = floor(timeinc/myh)
myk = myh^(0.55)
xylimit = 1.5*max(abs(mydata[,1:2]))
M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)
mm = length(xvec)

xcord = mydata[,1]
ycord = mydata[,2]

# create grid of theta values at which we'd like to evaluate likelihood
theta12grid = c(-5:5)*0.1 + 1
n12 = length(theta12grid)
theta3grid = c(-5:5)*0.2 + 2
n3 = length(theta3grid)

# pomp initialization
source('pompinit.R')
numparticles = 1000

# now evaluate likelihood using both DTQ and pomp
dtqlik = array(dim=c(n12,n12,n3))
pomplik = array(dim=c(n12,n12,n3))

for (ii in c(1:n12))
{
  for (jj in c(1:n12))
  {
    for (kk in c(1:n3))
    {
      thetavec = truethetavec
      thetavec[1] = theta12grid[ii]
      thetavec[2] = theta12grid[jj]
      thetavec[3] = theta3grid[kk]
      dtqtmp = Rdtq2d(thetavec, xcord, ycord, h=myh, numsteps=myns, k=myk, yM=xylimit)
      dtqtmp[dtqtmp <= 2.2e-16] = 2.2e-16
      dtqlik[ii,jj,kk] = sum(log(dtqtmp))

      myparams = c(theta1=thetavec[1],theta2=thetavec[2],theta3=thetavec[3],
                   X1.0 = mymod.dat[1,]$Y1, X2.0 = mymod.dat[1,]$Y2)
      pomppf = pfilter(mymod, Np = numparticles, params = myparams)
      pomplik[ii,jj,kk] = logLik(pomppf)
      print(c(thetavec,dtqlik[ii,jj,kk],pomplik[ii,jj,kk]))
      flush.console()
    }
  }
}
save.image(file = 'likcompare.RData')


