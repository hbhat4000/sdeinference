# This code generate bootstrap samples of inverted parameters theta and use the
# Hermite function expansion to generate samples for drift functions

library('boot')
source('optimtestfun.R')
source('myhermfun.R')

theta.boot <- boot(fd ,myoptfun, R=4)
bootsamp = theta.boot$t
d = dim(bootsamp)


xvec = seq(-5,5,0.0001)
npts = length(xvec)
numpara = d[2]-1
hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)

sum = 0
for (i in c(1:numpara))
{
	bootdriftmat = outer(bootsamp[,i], hermitefuncmat[,i])
	bootdriftmat = sum  + bootdriftmat
}