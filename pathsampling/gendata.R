# clear memory
rm(list=ls(all=TRUE))

# load required package
library('orthopolynom')

# steps required to get this code to work
# 1) create an artificial drift function and a trajectory
# 2) write main EM loop
# 2a) E-step: use diffusion bridges to get many fine paths
# 2b) M-step: using entire collection of paths, solve regression problem
# 3) rinse, lather, repeat

# this code ONLY implements 1) from above!

# use orthogonal polynomials
# set number of degrees of freedom (numdof = max degree + 1)
numdof = 4

# sample frozen set of coefficients ~ U(-1,1)
# set.seed(123)
# frozencoeff = 2*runif(n=numdof) - 1
# frozencoeff = frozencoeff/10
frozencoeff = c(1,2,-1,-.5)

# suppose diffusion coeff is known
g = 1/2

# create frozen drift function
normalized.p.list = hermite.he.polynomials(n=numdof-1, normalized=TRUE)
frozendrift <- function(x)
{
  polyvals = matrix(nrow=length(x),unlist(polynomial.values(normalized.p.list,x)))
  as.numeric(polyvals %*% frozencoeff)
}

# create trajectory
ic = 0
ft = 10.0
numsteps = 25000  # total number of "fine" steps to take
savesteps = 100  # total number of times to save the solution
traj = numeric(length=savesteps)
tvec = numeric(length=savesteps)
h = ft/numsteps
h12 = sqrt(h)
x = ic
j = 1
for (i in c(1:numsteps))
{
  x <- x + frozendrift(x)*h + g*h12*rnorm(n=1)
  if ((i %% (numsteps/savesteps)) == 0)
  {
    traj[j] = x
    tvec[j] = i*h
    j = j + 1
  }
}
traj = c(ic,traj)
tvec = c(0,tvec)

# save the data
save(tvec,traj,file='nem5.RData')


