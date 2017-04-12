rm(list = ls(all = TRUE))

library("matrixcalc")

# xtraj = matrix(nrow = 1, ncol = 2)

# load necessary functions
source('dtq_main.R')
source('Dtheta.R')


k = 0.01
M = 800
deltat = 1
numsteps = 50
h = deltat/numsteps
theta = c(1, 0, 2)
init = 1
final = 3

source('kolmogorov_compare.R')

if(numsteps >= 1)
{
  completelik_front = dtq_complete_front(theta, h, k, M, numsteps, init, final)
  completelik_back = dtq_complete_back(theta, h, k, M, numsteps, init, final)
  
  # compare against exact solution
  exactcompletelik = transition_forward(theta,x=final,y=init,t=deltat)
  print(c(completelik_front$lik, completelik_back,exactcompletelik))
}

if(numsteps >= 2)
{
  # first step should be equal to the last step if 2 steps
  firststeplik = dtq_firststep(theta, h, k, M, numsteps, init, final)
  ourfirststep = firststeplik / completelik_back
  
  # check normalization
  print(sum(ourfirststep)*k)
  
  part1 = transition_forward(theta,x=grid,y=init,t=h)
  part2 = transition_forward(theta,x=final,y=grid,t=(deltat-h))
  exactfirststeplik = log(hadamard.prod(part1,part2)/exactcompletelik)
  
  laststeplik = dtq_laststep(theta, h, k, M, numsteps, init, final)
  ourlaststep = laststeplik / completelik_back

  # check normalization
  print(sum(ourlaststep)*k)
  
  part1 = transition_forward(theta,x=grid,y=init,t=(deltat-h))
  part2 = transition_forward(theta,x=final,y=grid,t=h)
  exactlaststeplik = log(hadamard.prod(part1,part2)/exactcompletelik)
  
  # par(mfrow=c(1,2))
  # plot(ourfirststep)
  # lines(exactfirststeplik,col='red')
  # plot(ourlaststep)
  # lines(exactlaststeplik,col='red')
}
# 
if(numsteps >= 3)
{
  # numsteps = F+1
  for (j in c(1:(numsteps-2))) {
    internallik = dtq_internal(theta, h, k, M, numsteps, init, final, j)
    internallik = internallik/completelik_back
    print(c(j, sum(internallik)*k^2))
  }
}
  
  
# nsamples = 10
# thetamat = matrix(nrow = (nsamples + 1), ncol = 3)
# # define log prior
# logprior <- function(z)
# {
#   return(sum(dnorm(x = z, mean = 0, sd = 100, log = TRUE)))
# }
# thetamat[1,] = theta
# 
# propZ <- function(n)
# {
#   return(rnorm(n, sd = 0.025))
# }
# 
# loglik = numeric(length = (nsamples + 1))
# logpost = numeric(length = (nsamples + 1))
# ar = numeric(length = nsamples)
# 
# rawlik = dtq(thetamat[1,], h, k, M, 1,  init = xtraj[1,1], final = xtraj[1,2])
# rawlik[rawlik <= 2.2e-16] = 0
# loglik[1] = sum(log(rawlik))
# logpost[1] = loglik[1] + logprior(thetamat[1,])
# 
# # Metropolis algorithm
# for (i in c(1:nsamples))
# {
#   thetastar = thetamat[i,] + propZ(n=1)
#   rawlik = dtq(thetastar, h, k, M, 1, xtraj[1,1], final = xtraj[1,2])
#   rawlik[rawlik <= 2.2e-16] = 0
#   thisloglik = sum(log(rawlik))
#   thislogpost = thisloglik + logprior(thetastar)
#   ratio = exp(thislogpost - logpost[i])
#   if (ratio > runif(n=1))
#   {
#     thetamat[i+1,] = thetastar
#     loglik[i+1] = thisloglik
#     logpost[i+1] = thislogpost
#     ar[i] = 1
#   }
#   else
#   {
#     thetamat[i+1,] = thetamat[i,]
#     loglik[i+1] = loglik[i]
#     logpost[i+1] = logpost[i]
#   }
#   print(mean(ar[1:i]))
# }
