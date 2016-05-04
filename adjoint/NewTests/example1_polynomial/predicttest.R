
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


# dX_t = theta1 4(X_t - X_t^3) dt + dW_t
# This function compute the prediction paths for the SDE with inverted parameters

predict <- function(fdtest, thetaapprox)
{

	
	numpara = length(thetaapprox)
	diffcoeff = thetaapprox[numpara]
	h = 0.0001
	littlet = 1
	bigt = 50


	nsteps = ceiling(bigt/h)
	nsaves = ceiling(bigt/littlet)
	hilt = ceiling(littlet/h)
	stopifnot((nsteps == (nsaves*hilt)))

	ntrials = 100
	x0 = 0
	h12 = sqrt(h)
	xtrajnew = matrix(0,nrow=ntrials,ncol=(nsaves+1))
	xtrajnew[,1] = rep(x0,times=ntrials)

	for (i in c(1:nsaves))
	{
    	# print loop counter
    	#print(i)
    	#flush.console()

    	x = xtrajnew[,i]

    	for (j in c(1:hilt))
		fapprox = driftfun(thetaapprox, x) 
        	x = x + fapprox*h + h12*diffcoeff*rnorm(n=ntrials)

    	xtrajnew[,(i+1)] = x
	}
	tvec = seq(from=0,to=bigt,by=littlet)
	xtrajnew = rbind(tvec,xtrajnew)
	return(xtrajnew)

}
