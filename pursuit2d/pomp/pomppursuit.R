rm(list=ls(all=TRUE))
library('pomp')
library('magrittr')
library('dplyr')

load('~/Desktop/newfakedata_fullres.RData')
mydata = X[seq(from=1,to=2001,by=40),]
mymod.dat = data.frame(t=mydata[,3],Y1=as.numeric(mydata[,1]),Y2=as.numeric(mydata[,2]))

step.fun <- Csnippet("
 double dW1 = rnorm(0,sqrt(dt));
 double dW2 = rnorm(0,sqrt(dt));
 X1 += -theta1*theta1*X2*dt + (0.16/(2*M_PI))*theta1*theta1*dW1;
 X2 += 2*M_PI*X1*dt + (0.16)*dW2;
")

mymod <- pomp(data=mymod.dat,time="t",t0=0,
              rprocess=euler.sim(step.fun=step.fun,delta.t=0.01),
              statenames=c("X1","X2"),paramnames=c("theta1"))

rmeas <- Csnippet("
 Y1 = X1 + rnorm(0,1e-2);
 Y2 = X2 + rnorm(0,1e-2);
")

mymod <- pomp(mymod,rmeasure=rmeas,statenames=c("X1","X2"))

dmeas <- Csnippet("
 lik = dnorm(Y1,X1,1e-2,give_log) + dnorm(Y2,X2,1e-2,give_log);
")

mymod <- pomp(mymod,dmeasure=dmeas,statenames=c("X1","X2"))

mymod  %<>%
  pomp(dprior=Csnippet("lik = dnorm(theta1,0,100,1);"),paramnames=c("theta1"))

startvec = c(1,mymod.dat[1,]$Y1,mymod.dat[1,]$Y2)
sdvec = c(1)
names(startvec) = c("theta1","X1.0","X2.0")
names(sdvec) = c("theta1")

mymod %>% pmcmc(Nmcmc=200,Np=200,start=startvec,proposal=mvn.rw.adaptive(rw.sd=sdvec,scale.start=100,shape.start=100)) -> chain
# chain %<>% pmcmc(Nmcmc=10000,proposal=mvn.rw(covmat(chain)))


