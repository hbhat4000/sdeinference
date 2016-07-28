rm(list=ls(all=TRUE))

library('Sim.DiffProc')

# generate simulated data
f = expression(1*x*(4 - x^2))
g = expression(0.5)
ntrials = 100
h = 0.0001
littlet = 1
bigt = 25
nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
xtraj = matrix(0,nrow=ntrials,ncol=(nsaves+1))
xtraj[,1] = rnorm(n=ntrials)
for (i in c(1:ntrials))
{
  sim = snssde1d(drift=f,diffusion=g,x0=xtraj[i,1],M=1,N=nsteps,t0=0,T=bigt,Dt=h)
  myind = (seq(from=1,by=(nsteps/nsaves),to=length(sim$X)))[-1]
  xtraj[i,2:(nsaves+1)] = sim$X[myind]
  print(i)
  flush.console()
}

save(xtraj,file='doublewell.RData')
