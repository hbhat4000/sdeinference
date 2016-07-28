rm(list=ls(all=TRUE))

library('Sim.DiffProc')
load('doublewell.RData')

# use all 100 time series
xtrajts = ts(data=t(xtraj),start=0,end=26,deltat=1)

# only use 1 time series
# xtrajts = ts(data=xtraj[50,],start=0,end=26,deltat=1)

fx = expression(theta[1]*x*(theta[2] - x^2))
gx = expression(theta[3])

# the only methods that work are "ozaki" and "shoji"
fitmod = fitsde(data=xtrajts,drift=fx,diffusion=gx,start=list(theta1=8,theta2=8,theta3=2),pmle="ozaki")

# next steps
# put this inside an MCMC loop
# generate posteriors
# compare with ours

