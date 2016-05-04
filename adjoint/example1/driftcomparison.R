driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

source('myhermfun.R')

#thetaapprox = c( 0.038493, -5.318236, 0.096413, -11.118269, -0.014055 )
#thetaapprox = c( 0.250357, -8.470558, 1.445716, -23.680692, 2.729792, -18.510258, 1.816888, -9.997309 )
#thetaapprox = c( 3.207790, -8.968125, 31.134077, -20.407300, 102.480341, 5.536307, 162.358304, 43.157791, 127.916599, 51.188084, 41.002527, 18.821167 )

# parallel run 1
# thetaapprox = c( 3.153656, -8.067596, 30.577911, -14.964415, 100.554822, 19.962084, 159.185275, 62.993084, 125.337884, 65.144399, 40.160603, 22.829154 )

# parallel run 2
# thetaapprox = c( 0.041191, -9.143040, 0.093259, -24.995346, 0.100000, -14.625772, 0.047683, 5.844716, -0.100000, 18.963042, 0.100000, 7.920510 )

# small run
# thetaapprox = c( 0.064075, -0.593441, -0.100000 )

dx = 0.01
xvec = seq(-5,5,dx)
npts = length(xvec)
numpara = length(thetaapprox)
hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)
fapprox = driftfun(thetaapprox, xvec, hermitefuncmat) 
ftrue = 4*(xvec - xvec^3)
plot(xvec, fapprox, type='l', col='red')
lines(xvec, ftrue )

error = sum((ftrue-fapprox)^2)*dx
print(error)


#[1] 316686.2 error with h=0.1 and 5 Hermite functions with random initial coefficients
#[1] 305881.1 error with h=0.1 and 8 Hermite functions with random initial coefficients




