# clear memory
rm(list=ls(all=TRUE))

# load required packages
library('orthopolynom')
library('sde')
library('parallel')

# load data
load('~/nem.RData')
n = length(tvec)
dt = tvec[2]-tvec[1]
zoomfac = 100
h = dt/zoomfac
bign = (n-1)*zoomfac + 1

# steps required to get this code to work
# 1) create an artificial drift function and a trajectory
# 2) write main EM loop
# 2a) E-step: use diffusion bridges to get many fine paths
# 2b) M-step: using entire collection of paths, solve regression problem
# 3) rinse, lather, repeat

# this code implements 2-3) from above, in parallel!

# use orthogonal polynomials
# set number of degrees of freedom (numdof = max degree + 1)
numdof = 3
normalized.p.list = hermite.he.polynomials(n=numdof-1, normalized=TRUE)
mypoly <- function(x)
{
  matrix(nrow=length(x),unlist(polynomial.values(normalized.p.list,x)))
}

# suppose diffusion coeff is known
g = 1/8

# set up main EM loop
# set EM exit tolerance, initial guess for theta
mytol = 1e-3
theta = rnorm(n=numdof)
done = 0
numpaths = 10
sigmafn = expression(g)
sigmaxfn = expression(0.0)
numiter = 0

while (done==0)
{
  numiter = numiter + 1
  print(numiter)
  
  # from the current theta vector, construct the drift fn as an expression
  driftfn = normalized.p.list[[1]]*theta[1]
  for (i in c(2:numdof))
    driftfn <- driftfn + normalized.p.list[[i]]*theta[i]
  
  # convert the polynomial to a legit R expression
  driftfn = parse(text=as.character(driftfn))
  
  # need to do this across all time steps
  # create a large number of trajectories
  
  # need to make more sample paths
  mmat = matrix(0,nrow=numdof,ncol=numdof)
  rvec = numeric(length=numdof)
  
  myfun <- function(i)
  {
    as.numeric(DBridge(x=traj[i],y=traj[i+1],t0=0,T=dt,
                       delta=h,drift=driftfn,sigma=sigmafn,sigma.x=sigmaxfn))
  }
  
  for (s in c(1:numpaths))
  {
    myfunevals = mclapply(X=c(1:(n-1)),FUN=myfun)
    completedata = unlist(myfunevals)[-(zoomfac+1)*c(1:(n-2))]
    pp = mypoly(completedata[-bign])
    mmat = mmat + h*(t(pp) %*% pp)/numpaths
    rvec = rvec + as.numeric(t(diff(completedata)) %*% pp)/numpaths
  }
  
  # plot original data in red (circles) and completed data in black (lines)
  # plot(tvec,traj,col='red')
  # lines(seq(from=0,to=tvec[26],by=h),completedata)
  
  # M-step
  newtheta = solve(mmat, rvec)

  check = sum(abs(newtheta - theta))
  if (check < mytol)
  {
    print(check)
    print(theta)
    done = 1
  }
  theta = newtheta
  print(theta)
}

