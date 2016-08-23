# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

rm(list=ls(all=TRUE))
# Creating test data
trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

# create independent x-values 
# x = c(-15:15)
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)

# Stochastic Equation: y = ax + b + noise
# create dependent values according to y = ax + b + N(0, sd)
# 31 dimensional space
y <-  trueA * x + trueB + rnorm(n = sampleSize, mean = 0, sd = trueSd)

plot(x, y, main="Test Data")

# likelihood function
likelihood <- function(param)
{
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  # generates log likelihood using dnorm
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

# Example: plot the likelihood profile of the slope a
slopevalues <- function(x)
{
  return(likelihood(c(x, trueB, trueSd)))
}

slopelikelihoods <- lapply(seq(3, 7, by = 0.05), slopevalues )
plot(seq(3, 7, by= 0.05), slopelikelihoods, type = "l", xlab = "values of slope parameter a", ylab = "Log likelihood")

# Prior distribution
prior <- function(param)
{
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min = 0, max = 10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min = 0, max = 30, log = T)
  return(aprior + bprior + sdprior)
}

# Posterior distribution
posterior <- function(param)
{
  return (likelihood(param) + prior(param))
}

######## Metropolis algorithm ################
proposalfunction <- function(param)
{
  return(rnorm(3, mean = param, sd = c(0.1, 0.5, 0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations)
{
  chain = array(dim = c(iterations+1, 3))
  chain[1,] = startvalue
  for (i in 1:iterations)
  {
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab)
    {
      chain[i+1,] = proposal
    }
    else
    {
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(4,0,10)
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1 - mean(duplicated(chain[-(1:burnIn),]))

### Summary: #######################

par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1], nclass = 30, main = "Posterior of a", xlab = "True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = trueA, col = "red" )
hist(chain[-(1:burnIn),2], nclass = 30, main = "Posterior of b", xlab = "True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = trueB, col = "red" )
hist(chain[-(1:burnIn),3], nclass = 30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = trueSd, col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )

# for comparison:
# summary(lm(y~x))
