driftfun <- function(c0, y) 
{
  # computes drift as a matrix
  f = (c0[1])*(c0[2] - y)
  return(f)
}

difffun <- function(c0, y)
{
  # replicates diffusion as a long vector
  g = rep(c0[3], length(y))
  return(g)
}

difffun1 <- function(c0, y)
{
  # replicates diffusion term as a matrix
  # dd1 = dim(y)[1]
  # dd2 = dim(y)[2]
  outmat = matrix(c0[3],nrow=nrow(y),ncol=ncol(y))
  return(outmat)
}

source('integrandmat.R')
source('Dtheta.R')

# This function returns the objective function and the gradient for a given parameter vector c0
# The gradient is computed using the adjoint method
# This function is written to compute the negative loglikelihood and the gradient for an SDE with 
# constant diffussion function given by the last component of the input paraeter vector c0
# The drift function f is considered as an expansion of Hermite functions
# The parameter vector c0 contains the coefficient of the Hermite function expansion of the drift f
# except the last element of c0

gradobjfun <- function(c0, h, k, bigm, littlet, data)
{
  numsteps = ceiling(littlet/h)
  xvec = c((-bigm):bigm)*k
  npts = length(xvec)
  ndata = ncol(data)
  numpara = length(c0)-1

  # hermitefuncmat is a matrix of Hermite functions evaluated on the grid xvec, stored along the columns
  # Along the rows the degree of the Hermite function increases starting from zero th degree to length(c0)-1
  # hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)

  # A is the kernel matrix K of the QDT method without the factor of grid spacing k
  A = integrandmat(xvec, xvec, h, driftfun, difffun1, c0)
    
  # D is a list of number of length(c0) matrices. Each matrix is the derivalive of the kernel matrix K w.r.t. one parameter.
  D = Dtheta(xvec, xvec, h, driftfun, difffun1, c0)

  nc0 = length(c0)
  initder = list(NULL)

  grad = 0
  cost = 0

  for (j in c(1:(ndata-1)))
  {
    # pdfmat contains all the pdfs at each littlet timesteps
    pdfmat = matrix(0, nrow = npts, ncol = (numsteps-1))

    pdfmat[,1] = rowMeans(integrandmat(xvec, data[2:nrow(data),j], h, driftfun, difffun1, c0))
    pdf = as.vector(pdfmat[,1])

    # solve the state equation
    startstep = 2
    for (curstep in c(startstep:(numsteps-1)))
    {
      pdf = k*A %*% pdf
      pdfmat[,curstep]=pdf
    }
    gammamat = k*integrandmat(data[2:nrow(data), (j+1)], xvec, h, driftfun, difffun1, c0)
    likelihood = gammamat %*% pdf

    # cost is the objective function given by the log likelihood without the negative sign
    cost = cost + sum(log(likelihood))

    # this is the last adjoint vector (u_j+(F-1)/F = y equation 23b  in notes)
    adj = apply(gammamat/as.vector(likelihood), 2, sum)

    # stores all the adjoint vectors u_j+(F-N)/F along columns for each internal timestep
    adjmat = matrix(0, nrow = (numsteps-1), ncol = npts) 
    adjmat[(numsteps-1),]=adj

    # solve the adjoint equation backward
    tA = t(A)
    for (curadjstep in c((numsteps-2):1))
    {
      adj = k*tA %*% adj
      adjmat[curadjstep,] = adj	
    }

    gradvec = 0*c(1:nc0)
    termvec = 0*c(1:nc0)

    initderivs = Dtheta(xvec, data[2:nrow(data),j], h, driftfun, difffun1, c0)
    DW = Dtheta(data[2:nrow(data),(j+1)], xvec, h, driftfun, difffun1, c0)

    # compute the two terms of the gradient equation (equation 24 in notes) w.r.t. each parameter
    for (i in c(1:nc0))
    {
      termvec[i] = -sum((k*DW[[i]] %*% pdf)/likelihood)
      M = -k*D[[i]] %*% pdfmat
      M[,2:(numsteps-1)] = M[,1:(numsteps-2)]
      # rowMeans(initderivs[[i]]) is the derivative of \hat(p)_j+1/F w.r.t i th parameter 
      M[,1] = -rowMeans(initderivs[[i]])
      gradvec[i] = sum(t(M)*adjmat)
    }

    # grad is the gradient of the objective function
    grad = grad + termvec + gradvec		
  }
  return(list(obj=-cost,grad=grad))
}


