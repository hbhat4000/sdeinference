
driftfun <- function(c0, y) 
{
  # computes drift as a matrix
  nc0 = length(c0) - 1
  if (class(y)=="numeric") y = as.matrix(y) 
  f = matrix(data = c0[1],nrow=nrow(y),ncol=ncol(y))
  for (i in c(2:nc0))
    f = f + c0[i]*y^(i-1)

  # f = c0[1]*y + c0[2]*y^3
  return(f)
}


source('predicttest.R')

compare <- function(thetaapprox, fdtest, fdtrain)
{
	dx = 0.00001
	xvec = seq(-5,5,dx)
	npts = length(xvec)
	numpara = length(thetaapprox)
	fapprox = driftfun(thetaapprox, xvec) 
	ftrue = 4*(xvec - xvec^3)
	plot(xvec, fapprox, type='l', col='red')
	lines(xvec, ftrue )

	error = sum((ftrue-fapprox)^2)*dx
	print(error)
	

	fdnew = predict(fdtest, thetaapprox)	
	dimvec = dim(fdnew)
	vec1 = c(1:ceiling(dimvec[2]/2))
	fdtest_approx = fdnew[-1,-vec1]
	mse_mat = (fdtest-fdtest_approx)^2	
	mse = apply(mse_mat,2,mean)
	return(mse)
}

	