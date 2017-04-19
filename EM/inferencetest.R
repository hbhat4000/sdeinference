rm(list = ls(all = TRUE))

library("matrixcalc")

# load necessary functions
source('dtq_main.R')
source('Dtheta.R')

# tvec, xtraj
load('fakedata.RData')
xtraj = xtraj[1:2, 1:10]
xlen = length(xtraj[2,])

k = 0.05
M = 160
deltat = 1
h = 0.1
numsteps = deltat / h
# numsteps = 20
# h = deltat/numsteps
theta = c(2, 1, 0.5)
theta0 = c(1, 1, 1)
# init = 1
# final = 3

# allout = dtq_all(theta, h, k, M, numsteps, init, final)

# compute log G(x,y,theta) and its gradient w.r.t. theta
logG <- function(x,y,theta,f,g,h)
{
  # first compute log G by itself
  evalf = f(theta, y, returngrad=TRUE)
  evalg = g(theta, y, returngrad=TRUE)
  fval = evalf$val
  gval = evalg$val
  part1 = (x - y - fval*h)
  lg = -0.5*log(2*pi*(gval^2)*h) - (part1^2)/(2*(gval^2)*h)
  
  # next compute d(log G)/df and d(log G)/dg
  dlgdf = part1/(gval^2)
  dlgdg = -1/gval + (part1^2)/((gval^3)*h)
  
  # now loop over the entries of theta and compute the gradient 
  n = length(theta)
  
  # two cases for gradlg
  if (class(x)=="numeric")
  {
    # case 1: only one of {x,y} is a vector and the other is a scalar
    gradlg = matrix(0,nrow=max(length(x),length(y)),ncol=n)
    for (i in c(1:n))
    {
      gradlg[,i] = dlgdf * evalf$grad[,i] + dlgdg * evalg$grad[i]
    }
  }
  else if (class(x)=="matrix")
  {
    # case 2: both x and y are "gridmat" matrices of the same dimension
    gradlg = array(0,dim=c(dim(x),n))
    for (i in c(1:n))
    {
      gradlg[,,i] = dlgdf * evalf$grad[,,i] + dlgdg * evalg$grad[i]
    }
  }
  return(list(lg=lg,gradlg=gradlg))
}

# goal: compute the Q function and its gradient w.r.t. theta
qfun <- function(theta, allout, f, g, h, k, M, numsteps, x)
{
  grid = c((-M):M)*k
  gridmat = replicate(length(grid), grid)
  init = x[1]
  final = x[2]
  q = 0
  n = length(theta)
  gradq = numeric(length=n)
  
  # first term in the summation (j=1 case)
  part1 = logG(x=grid,y=init,theta,f,g,h)
  q = q + sum(part1$lg * allout$first)*k
  
  # loop over components of theta
  for (i in c(1:n))
    gradq[i] = gradq[i] + sum(part1$gradlg[,i] * allout$first)*k
  
  # last term in the summation (j=F case)
  part2 = logG(x=final,y=grid,theta,f,g,h)
  q = q + sum(part2$lg * allout$last)*k
  
  # loop over components of theta
  for (i in c(1:n))
    gradq[i] = gradq[i] + sum(part2$gradlg[,i] * allout$last)*k
  
  # all intermediate terms
  part3 = logG(x=gridmat,y=t(gridmat),theta,f,g,h)
  for (j in c(1:(numsteps-2)))
  {
    q = q + sum(part3$lg * allout$pdf2d[[j]])*k^2
    for (i in c(1:n))
      gradq[i] = gradq[i] + sum(part3$gradlg[,,i] * allout$pdf2d[[j]])*k^2
  }
  return(list(objective=-q,gradient=-gradq))
}

simpqfun <- function(theta, thistheta)
{
  theta0 = thistheta
  obj = 0
  grad = 0
  for (i in c(1:(xlen - 1)))
  {
    init = xtraj[2,i]
    final = xtraj[2, i+1]
    
    allout = dtq_all(theta0, h, k, M, numsteps, init, final)
    qinterval = qfun(theta, allout=allout, f, g, h, k, M, numsteps, x=c(init,final))
    obj = obj + qinterval$objective
    grad = grad + qinterval$gradient
    qinterval = list(objective = obj, gradient = grad)
  }
  return (qinterval)
}

library('nloptr')

mylb = c(0, 0, 0.01)
myub = c(10, 10, 2)

# Q: what value of \theta should each outer maximization start from?
# currently, \theta parameter always gets initialized to c(2, 1, 0.5)
for(i in c(1:5)) 
{
  print(c(theta, theta0))
  # theta = theta0
  res = nloptr(x0 = theta, eval_f = simpqfun, opts = list("algorithm"="NLOPT_LD_LBFGS", "print_level" = 3, "check_derivatives"=TRUE, "xtol_abs" = 1e-8), lb = mylb, ub = myub, thistheta = theta0)
  
  # \theta^{k+1} = \arg \max_{\theta} Q(\theta | \theta^{k})
  theta0 = res$solution
}

