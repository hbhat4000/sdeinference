library('pomp')
library('magrittr')
library('dplyr')

mymod.dat = data.frame(t=mydata[,3],Y1=as.numeric(mydata[,1]),Y2=as.numeric(mydata[,2]))

step.fun <- Csnippet("
 double dW1 = rnorm(0,sqrt(dt));
 double dW2 = rnorm(0,sqrt(dt));
 X1 += theta1*theta1*X2*dt + (0.1+exp(-X1*X1)*9/16)*dW1;
 X2 += -theta2*theta2*X1*dt + theta3*theta3*(1-X1*X1)*X2*dt + (0.1+exp(-X2*X2)*9/16)*dW2;
")

mymod <- pomp(data=mymod.dat,time="t",t0=0,
              rprocess=euler.sim(step.fun=step.fun,delta.t=myh),
              statenames=c("X1","X2"),paramnames=c("theta1","theta2","theta3"))

rmeas <- Csnippet("
 Y1 = X1 + rnorm(0,1e-2);
 Y2 = X2 + rnorm(0,1e-2);
")

mymod <- pomp(mymod,rmeasure=rmeas,statenames=c("X1","X2"))

dmeas <- Csnippet("
  lik = dnorm(Y1,X1,1e-2,1) + dnorm(Y2,X2,1e-2,1);
")

mymod <- pomp(mymod,dmeasure=dmeas,statenames=c("X1","X2"))
