rm(list=ls(all=TRUE))

library('pomp')
library('magrittr')
library('dplyr')
library('foreach')
library('doParallel')
registerDoParallel(cores=24)

tvec = c(0:25)
load('doublewell.RData')
xvec = xtraj[1,]
mymod.dat = data.frame(t=tvec,X=as.numeric(xvec))

step.fun <- Csnippet("
  double dW = rnorm(0,sqrt(dt));
  X += theta1*X*(theta2-X*X)*dt + theta3*dW;
")

mymod <- pomp(data=mymod.dat,time="t",t0=0,
              rprocess=euler.sim(step.fun=step.fun,delta.t=0.01),
              statenames="X",paramnames=c("theta1","theta2","theta3"))

# simStates = simulate(mymod,nsim=10,params=c(theta1=.5,theta2=1,theta3=.25,X.0=0),states=TRUE)

# by itself, all the following does is compute the log likelihood
pf <- pfilter(mymod,Np=10000,params=c(theta1=1.0,theta2=0.1,theta3=0.25,X.0=0),save.states=TRUE)

# mymod  %<>% 
#   pomp(dprior=Csnippet("
#     lik = dnorm(theta1,0.5,1,1) + dnorm(theta2,2.0,10,1) + dexp(exp(sigeps),1,1);
#     lik = (give_log) ? lik : exp(lik);
#   "), paramnames=c("theta1","theta2","sigeps"))

# startvec = c(1,0.1,0,as.numeric(xvec[1]))
# sdvec = c(0.05,0.05,0.02)
# names(startvec) = c("theta1","theta2","sigeps","X.0")
# names(sdvec) = names(startvec)[1:3]
# starts = data.frame(matrix(rep(startvec,each=200),ncol=3))
# colnames(starts) = c("theta1","theta2","sigeps")

# mymod %>% pmcmc(Nmcmc=2000,Np=200,start=startvec,proposal=mvn.rw.adaptive(rw.sd=sdvec,scale.start=100,shape.start=100)) -> chain
# chain %<>% pmcmc(Nmcmc=10000,proposal=mvn.rw(covmat(chain)))

# save(chain,file="chains_x0.RData")

