rm(list = ls(all = TRUE))

source('dtq_with_grad.R')

thetavec = c(0,0.8,0,-0.2,0.4)

h = 0.0001
littlet = 1
bigt = 25

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

ntrials = 300
h12 = sqrt(h)
xtraj = matrix(0, nrow = ntrials, ncol = (nsaves + 1))

# initial condition centered at origin so both -ve and +ve initial conditions generated
xtraj[,1] = rnorm(n = ntrials, mean = 0, sd = 1) 

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x = xtraj[,i]

    for (j in c(1:hilt))
        x = x + driftfun(thetavec,x)*(h) + (h12)*difffun(thetavec,x)*(rnorm(n = ntrials))      

    xtraj[,(i+1)] = x
}

browser()

if (sum(is.na(xtraj[,(nsaves+1)])) > 0)
    xtraj = xtraj[-which(is.na(xtraj[,(nsaves+1)])),]

tvec = seq(from = 0, to = bigt, by = littlet)
xtraj = rbind(tvec, xtraj)
save(xtraj, file = 'fakedata.RData')


