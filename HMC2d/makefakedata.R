rm(list = ls(all = TRUE))

source('driftdiff.R')
source('truethetavec.R')

h = 0.0001
t = 1
T = 25

nsteps = ceiling(T/h)   # 250000
nsaves = ceiling(T/t)   # 25
hilt = ceiling(t/h)		# 1000
stopifnot((nsteps == (nsaves*hilt)))

theta = truethetavec
ntrials = 1
h12 = sqrt(h)

xtraj1 = matrix(0, nrow = ntrials, ncol = (nsaves + 1))
xtraj1[,1] = rnorm(n = ntrials, mean = 0, sd = sqrt(h)) 
xtraj2 = matrix(0, nrow = ntrials, ncol = (nsaves + 1))
xtraj2[,1] = rnorm(n = ntrials, mean = 0, sd = sqrt(h)) 

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    x1 = xtraj1[,i]
    x2 = xtraj2[,i]
    x = c(x1, x2)
    
    for (j in c(1:hilt))
    {
        # x1 = x1 + (f1(x1, x2, theta))*(h) + (h12)*(g1(x1, x2, theta))*(rnorm(n = ntrials))      
        # x2 = x2 + (f2(x1, x2, theta))*(h) + (h12)*(g2(x1, x2, theta))*(rnorm(n = ntrials))
        x1 = x1 + (f1(x, theta))*(h) + (h12)*(g1(x, theta))*(rnorm(n = ntrials))      
        x2 = x2 + (f2(x, theta))*(h) + (h12)*(g2(x, theta))*(rnorm(n = ntrials))
    }

    xtraj1[,(i+1)] = x1
    xtraj2[,(i+1)] = x2
}

tvec = seq(from = 0, to = T, by = t)
xtraj = list(tvec, xtraj1, xtraj2)
save(xtraj, file = 'fakedata.RData')