rm(list = ls(all = TRUE))

# solve inverse problem for linear SDE
# dX_t = theta1 (theta2 - X_t) dt +  theta3 dW_t
thetavec = c(0.5, 1, 0.25)

h = 1e-6
littlet = 2e-1
bigt = 500

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = 2e5 # ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

ntrials = 1
h12 = sqrt(h)
xtraj = matrix(0, nrow = ntrials, ncol = (nsaves + 1))

# initial condition centered at origin so both -ve and +ve initial conditions generated
xtraj[,1] = rnorm(n = ntrials) 

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x = xtraj[,i]

    for (j in c(1:hilt))
        x = x + (thetavec[1]*(thetavec[2] - x))*(h) + (h12)*(thetavec[3])*(rnorm(n = ntrials))

    xtraj[,(i+1)] = x
}

xtraj = xtraj + matrix(rnorm(n=((nsaves+1)*ntrials),mean=0,sd=0.1),nrow=ntrials,ncol=(nsaves+1))

xvec = as.numeric(xtraj)
tvec = seq(from=0,to=bigt,by=littlet)
write.table(matrix(tvec,nrow=1),file="tvec.csv",row.names=FALSE,col.names=FALSE,sep=',')
write.table(matrix(xvec,nrow=1),file="xvec.csv",row.names=FALSE,col.names=FALSE,sep=',')

# xtraj = rbind(tvec, xtraj)
# save(xtraj, file = 'fakedata1.RData')
