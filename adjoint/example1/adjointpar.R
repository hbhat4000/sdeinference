source('integrandmat.R')
source('Dtheta.R')
source('myhermfun.R')

library('foreach')
library('doParallel')
registerDoParallel(cores=24)

driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

difffun <- function(diffcoeff, y)
{
  # replicates diffusion term as a matrix
  dd1 = dim(y)[1]
  dd2 = dim(y)[2]
  vec = rep(diffcoeff, dd1)
  return(replicate(dd2, vec)) 
}

# This function returns the objective function and the gradient for a given parameter vector c0
# The gradient is computed using the adjoint method
# This function is written to compute the negative loglikelihood and the gradient for an SDE with constant diffussion function given by the input diffcoeff
# The drift function f is considered as an expansion of Hermite functions
# The parameter vector c0 contains the coefficient of the Hermite function expansion of the drift f

gradobjfun <- function(c0, h, k, bigm, littlet, data, diffcoeff)
{
    numsteps = ceiling(littlet/h)
    xvec = c((-bigm):bigm)*k
    npts = length(xvec)
    ndata = ncol(data)
    numpara = length(c0)

    # hermitefuncmat is a matrix of Hermite functions evaluated on the grid xvec, stored along the columns
    # Along the rows the degree of the Hermite function increases starting from zero th degree to length(c0)-1
    hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)

    # A is the kernal matrix K of the QDT method without the factor of grid spacing k
    A = integrandmat(xvec, xvec, h, driftfun, difffun, c0, hermitefuncmat, diffcoeff)
    
    # D is a list of number of length(c0) matrices. Each matrix is the derivalive of the kernal matrix K w.r.t. one parameter.
    D = Dtheta(xvec, xvec, h, driftfun, difffun, c0,  hermitefuncmat, diffcoeff)

    nc0 = length(c0)
    initder = list(NULL)

    # grad = 0
    # cost = 0

    allresults = foreach (j=(1:(ndata-1))) %dopar%
    {
		# pdfmat contains all the pdfs at each littlet timesteps
		pdfmat = matrix(0, nrow = npts, ncol = (numsteps-1))

		
		hermfuninit = myhermfun(numterms = numpara, xgrid = data[2:nrow(data),j])
		pdfmat[,1] = rowMeans(integrandmat(xvec, data[2:nrow(data),j], h, driftfun, difffun, c0, hermfuninit, diffcoeff))
		pdf = as.vector(pdfmat[,1])

		# solve the state equation
		startstep = 2
    		for (curstep in c(startstep:(numsteps-1)))
    		{
        		pdf = k*A %*% pdf
			pdfmat[,curstep]=pdf
        		
    		}
		
		gammamat = k*integrandmat(data[2:nrow(data), (j+1)], xvec, h, driftfun, difffun, c0, hermitefuncmat, diffcoeff)
            likelihood = gammamat %*% pdf

		# cost is the objective function given by the log likelihood without the negative sign
		cost = sum(log(likelihood))

		adj = apply(gammamat/as.vector(likelihood), 2, sum) # this is the last adjoint vector (u_j+(F-1)/F = y equation 23b  in notes)
            adjmat = matrix(0, nrow = (numsteps-1), ncol = npts) # stores all the adjoint vectors u_j+(F-N)/F along columns for each internal timestep
		adjmat[(numsteps-1),]=adj

		# solve the adjoint equation backward
    		for (curadjstep in c((numsteps-2):1))
    		{
        		adj = k*t(A) %*% adj
			adjmat[curadjstep,] = adj	
    		}

		gradvec = 0*c(1:nc0)
		termvec = 0*c(1:nc0)
		
		initderivs = Dtheta(xvec, data[2:nrow(data),j], h, driftfun, difffun, c0, hermfuninit, diffcoeff)

		DW = Dtheta(data[2:nrow(data),(j+1)], xvec, h, driftfun, difffun, c0, hermitefuncmat, diffcoeff)
		
		# compute the two terms of the gradient equation (equation 24 in notes) w.r.t. each parameter
		for (i in c(1:nc0))
		{
			termvec[i] = -sum((k*DW[[i]] %*% pdf)/likelihood)
			M = -k*D[[i]] %*% pdfmat
			M[,2:(numsteps-1)] = M[,1:(numsteps-2)]
            	M[,1] = -rowMeans(initderivs[[i]])# rowMeans(initderivs[[i]]) is the derivative of \hat(p)_j+1/F w.r.t i th parameter 
			gradvec[i] = sum(t(M)*adjmat)

		}

		# grad is the gradient of the objective function
		grad = termvec + gradvec		
                list(c=cost,g=grad)
    }
    cost = 0
    grad = 0
    for (i in c(1:length(allresults)))
    {
        cost = cost + allresults[[i]]$c
        grad = grad + allresults[[i]]$g
    }

    return(list(obj=-cost,grad=grad))
}


