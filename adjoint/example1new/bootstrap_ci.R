# This code generate bootstrap samples of inverted parameters theta and use the Hermite function expansion to generate samples for drift functions
library('boot')
source('optimtestfun.R')
source('myhermfun.R')

driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}


theta.boot <- boot(fd ,myoptfun, R=4)
bootsamp = theta.boot$t
d = dim(bootsamp)


xvec = seq(-5,5,0.0001)
npts = length(xvec)
numpara = d[2]-1
hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)

bootdrift = matrix(0, nrow = d[1], ncol = npts)
for (i in c(1:d[1]))
{
	c0 = bootsamp[i, (1:numpara)]
	bootdrift[i, ] = driftfun(c0, xvec, hermitefuncmat)
}







