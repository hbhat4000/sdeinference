rm(list = ls(all = TRUE))

# solve inverse problem for linear SDE
# dX_t = theta1 (theta2 - X_t) dt +  theta3 dW_t
thetavec = c(0.5, 1, 0.25)

h = 1e-6
littlet = 2e-1
bigt = 50

nsaves = ceiling(bigt/littlet)

ntrials = 1
h12 = sqrt(h)

# initial condition centered at origin so both -ve and +ve initial conditions generated
xvec = numeric(length = (nsaves+1))
tvec = numeric(length = (nsaves+1))
xvec[1] = rnorm(n=1)
tvec[1] = 0

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x = xvec[i]

    hilt = sample(c(4e4:(4e5-4e4)),size=1)

    for (j in c(1:hilt))
        x = x + (thetavec[1]*(thetavec[2] - x))*(h) + (h12)*(thetavec[3])*(rnorm(n = ntrials))

    tvec[i+1] = tvec[i] + hilt*h
    xvec[i+1] = x
}

xvec = xvec + rnorm(n=(nsaves+1),mean=0,sd=0.1)

write.table(matrix(tvec,nrow=1),file="tvec.csv",row.names=FALSE,col.names=FALSE,sep=',')
write.table(matrix(xvec,nrow=1),file="xvec.csv",row.names=FALSE,col.names=FALSE,sep=',')

# xvec = rbind(tvec, xvec)
# save(xvec, file = 'fakedata1.RData')
