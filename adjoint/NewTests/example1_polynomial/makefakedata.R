rm(list=ls(all=TRUE))

# solve inverse problem for nonlinear SDE
# dX_t = theta1 4(X_t - X_t^3) dt + dW_t


diffcoeff = 1.0

h = 0.00001
littlet = 0.5
bigt = 50

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

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
        x = x + 4*(x - x^3)*h + h12*diffcoeff*rnorm(n=ntrials)

    xtraj[,(i+1)] = x
}

tvec = seq(from=0,to=bigt,by=littlet)
xtraj = rbind(tvec,xtraj)
save(xtraj,file='fakedata_tpoint5_T50_nt100.RData')

