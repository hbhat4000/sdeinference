driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

source('myhermfun.R')
source('predicttest.R')

compare <- function(thetaapprox, fdtest, fdtrain, diffcoeff)
{

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
	
	fdnew = predict(fdtest, thetaapprox, diffcoeff)	
	dimvec = dim(fdnew)
	vec1 = c(1:ceiling(dimvec[2]/2))
	fdtest_approx = fdnew[-1,-vec1]
	mse_mat = (fdtest-fdtest_approx)^2	
	mse = apply(mse_mat,2,mean)
	return(mse)
}




