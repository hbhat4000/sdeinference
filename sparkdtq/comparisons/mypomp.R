library('pomp')

tvec = read.csv('../actual_bigrun12/tvec.csv',header=FALSE)
xvec = read.csv('../actual_bigrun12/xvec.csv',header=FALSE)
mymod.dat = data.frame(t=as.numeric(tvec),Y=as.numeric(xvec))
# mymod.dat = mymod.dat[1:21,]

step.fun <- Csnippet("
  double dW = rnorm(0,sqrt(dt));
  X += theta1*(theta2-X)*dt + theta3*dW;
")

mymod <- pomp(data=mymod.dat,time="t",t0=0,
              rprocess=euler.sim(step.fun=step.fun,delta.t=0.01),
              statenames="X",paramnames=c("theta1","theta2","theta3"))

# simStates = simulate(mymod,nsim=10,params=c(theta1=.5,theta2=1,theta3=.25,X.0=0),states=TRUE)

rmeas <- Csnippet("
  Y = X + rnorm(0,sigeps);
")

mymod <- pomp(mymod,rmeasure=rmeas,statenames="X",paramnames="sigeps")

dmeas <- Csnippet("
  lik = dnorm(Y,X,sigeps,give_log=1);
")

mymod <- pomp(mymod,dmeasure=dmeas,statenames="X",paramnames="sigeps")

pf <- pfilter(mymod,Np=10000,params=c(theta1=.5,theta2=1,theta3=.25,sigeps=0.1,X.0=-1.350519),save.states=TRUE)

# --> to do:
# figure out how to do perform MCMC inference
# try both particle and ABC methods
# compare inference for both parameters and states

