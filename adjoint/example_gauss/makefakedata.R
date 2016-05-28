rm(list=ls(all=TRUE))

# solve inverse problem for nonlinear SDE
# dX_t = theta1 exp((X_t-mu)^2/2 sigma^2)/sqrt(2 pi sigma^2) dt + dW_t


diffcoeff = 1.0

h = 0.0001
littlet = 0.5
bigt = 25

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

mu = 0
sigma = 1
ntrials = 100
x0 = 0
h12 = sqrt(h)
xtraj = matrix(0,nrow=ntrials,ncol=(nsaves+1))
xtraj[,1] = rep(x0,times=ntrials)

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x = xtraj[,i]

    for (j in c(1:hilt))
    {
	  zz = rnorm(n=ntrials)
        x = x + (exp(-(x-mu)^2/(2*sigma^2))/sqrt(2*pi*sigma^2))*h + h12*diffcoeff*zz
    }
    xtraj[,(i+1)] = x
}

tvec = seq(from=0,to=bigt,by=littlet)
xtraj = rbind(tvec,xtraj)
save(xtraj,file='fakedata_pointfive_100.RData')

