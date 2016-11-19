rm(list=ls(all=TRUE))
library('pomp')
library('magrittr')
library('dplyr')

# for reproducibility
set.seed(1)

load('newfakedata_fullres.RData')
mydata = X[seq(from=1,to=2001,by=20),]
mymod.dat = data.frame(t=mydata[,3],Y1=as.numeric(mydata[,1]),Y2=as.numeric(mydata[,2]))

ptm = proc.time()

timeinc = mydata[2,3] - mydata[1,3]
myh = timeinc/4

step.fun <- Csnippet("
 double dW1 = rnorm(0,sqrt(dt));
 double dW2 = rnorm(0,sqrt(dt));
 X1 += -theta1*theta1*X2*dt + (0.16/(2*M_PI))*theta1*theta1*dW1;
 X2 += theta2*theta2*X1*dt + (0.16/(2*M_PI))*theta2*theta2*dW2;
")

mymod <- pomp(data=mymod.dat,time="t",t0=0,
              rprocess=euler.sim(step.fun=step.fun,delta.t=myh),
              statenames=c("X1","X2"),paramnames=c("theta1","theta2"))

rmeas <- Csnippet("
 Y1 = X1 + rnorm(0,1e-2);
 Y2 = X2 + rnorm(0,1e-2);
")

mymod <- pomp(mymod,rmeasure=rmeas,statenames=c("X1","X2"))

dmeas <- Csnippet("
  lik = dnorm(Y1,X1,1e-2,1) + dnorm(Y2,X2,1e-2,1);
")

mymod <- pomp(mymod,dmeasure=dmeas,statenames=c("X1","X2"))

mymod  %<>%
  pomp(dprior=Csnippet("
    lik = dnorm(theta1,0,100,1) + dnorm(theta2,0,100,1);
  "),paramnames=c("theta1","theta2"))

startvec = c(1,1,mymod.dat[1,]$Y1,mymod.dat[1,]$Y2)
sdvec = c(1,1)
names(startvec) = c("theta1","theta2","X1.0","X2.0")
names(sdvec) = c("theta1","theta2")

numsamp=1000
burnin=100
nparticles=5
# rw.var = matrix(c(1,1,1,1), nrow=2, ncol = 2)

mymod %>% pmcmc(Nmcmc=200,Np=nparticles,start=startvec,proposal=mvn.rw.adaptive(rw.sd=sdvec,scale.start=.5,shape.start=1, target=0.234)) -> chain
chain %<>% pmcmc(Nmcmc=(numsamp+burnin),proposal=mvn.rw(covmat(chain)))
# chain %<>% pmcmc(Nmcmc=(numsamp+burnin),proposal=mvn.rw(rw.var))

# we will get exactly numsamp samples
mysamp = as.matrix(conv.rec(chain,"theta1"))[(burnin+2):(numsamp+burnin+1),1]

# this is because theta1^2 = 1/L, which should be 2*pi
print(summary(mysamp^2))

print(proc.time() - ptm)

