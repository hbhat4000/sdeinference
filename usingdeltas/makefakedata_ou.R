rm(list=ls(all=TRUE))

# solve inverse problem for nonlinear SDE
# dX_t = theta1(theta2- X_t) dt +  theta3 dW_t

theta1 = 0.5
theta2 = 1.5
theta3 = 1.0

h = 0.0001
littlet = 1
bigt = 25

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

ntrials = 500
h12 = sqrt(h)
xtraj = matrix(0,nrow=ntrials,ncol=(nsaves+1))
xtraj[,1] = rnorm(n=ntrials, mean=0, sd=0.2)

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x = xtraj[,i]

    for (j in c(1:hilt))
        x = x + (theta1)*(theta2 - x)*h + h12*theta3*rnorm(n=ntrials)

    xtraj[,(i+1)] = x
}

tvec = seq(from=0,to=bigt,by=littlet)
xtraj = rbind(tvec,xtraj)
save(xtraj,file='fakedata_ou_noise.RData')

