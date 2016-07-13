source('myhermfun.R')
# dX_t = theta1 4(X_t - X_t^3) dt + dW_t

driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

predict <- function(fdtest, thetaapprox, diffcoeff)
{
	numpara = length(thetaapprox)
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
	xtrajnew = matrix(0, nrow = ntrials, ncol = (nsaves+1))
	xtrajnew[,1] = rep(x0, times = ntrials)

	for (i in c(1:nsaves))
	{
    	x = xtrajnew[,i]

    	for (j in c(1:hilt))
    	{
    		hermitefuncmat = myhermfun(numterms = numpara, xgrid = x)
			fapprox = driftfun(thetaapprox, x, hermitefuncmat) 
        	x = x + fapprox*h + h12*diffcoeff*rnorm(n = ntrials)
        }

    	xtrajnew[,(i+1)] = x
	}

	tvec = seq(from = 0, to = bigt, by = littlet)
	xtrajnew = rbind(tvec, xtrajnew)
	return(xtrajnew)
}
