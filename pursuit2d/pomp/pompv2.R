rm(list=ls(all=TRUE))
library('pomp')
library('magrittr')
library('dplyr')

# for reproducibility
set.seed(1)

load('newfakedata_fullres.RData')
mydata = X[seq(from=1,to=2001,by=40),] # data in matrix form
mymod.dat = data.frame(Y1=as.numeric(mydata[,1]),Y2=as.numeric(mydata[,2]), t=mydata[,3])

step.fun <- Csnippet("
 double dW1 = rnorm(0,sqrt(dt));
 double dW2 = rnorm(0,sqrt(dt));
 X1 += -theta1*theta1*X2*dt + (0.16/(2*M_PI))*theta1*theta1*dW1;
 X2 += theta2*theta2*X1*dt + (0.16/(2*M_PI))*theta2*theta2*dW2;
")

mymod <- pomp(data=mymod.dat,time="t",t0=0,
              rprocess=euler.sim(step.fun=step.fun,delta.t=0.05),
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

burnin = 100
numsteps = 1000
totsteps = numsteps + burnin

artrack = numeric(length=(totsteps-1))
x = matrix(0, nrow = totsteps, ncol = 2)
x[1,1] = 2
x[1,2] = 2

oldpf <- pfilter(mymod, Np = 1000, params = c(theta1 = x[1,1], theta2 = x[1,2], X1.0 = mymod.dat[1,]$Y1, X2.0 = mymod.dat[1,]$Y2))
oldden <- logLik(oldpf)
print(oldden)

for(i in c(1:(totsteps-1))) {
	z = rnorm(n = 2, sd = 0.1)
	prop = x[i,] + z

	proppf <- pfilter(mymod, Np = 1000, params = c(theta1 = prop[1], theta2 = prop[2], X1.0 = mymod.dat[1,]$Y1, X2.0 = mymod.dat[1,]$Y2))
	propden <- logLik(proppf)
	print(propden)

	rho = exp(propden - oldden)

	u = runif(n = 1)
	if(rho > u) {
		x[i+1,] = prop
		oldpf = proppf
		oldden = propden
		artrack[i] = 1
	}
	else {
		x[i+1,] = x[i,]
		artrack[i] = 0
	}
}