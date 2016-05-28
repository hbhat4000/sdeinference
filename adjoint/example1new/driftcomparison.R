driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

source('myhermfun.R')
thetaapprox = c( -0.132798, -2.069660, -0.665487, -4.189549, -0.653228, 0.671999 )

xvec = seq(-5,5,0.0001)
npts = length(xvec)
numpara = length(thetaapprox)-1
hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)
fapprox = driftfun(c0=thetaapprox[1:numpara], xvec, hermitefuncmat) 
ftrue = 4*(xvec - xvec^3)
plot(xvec, fapprox, type='l', col='red')
lines(xvec, ftrue )