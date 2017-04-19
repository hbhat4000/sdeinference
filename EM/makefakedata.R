rm(list = ls(all = TRUE))
set.seed(1)
thetavec = c(1, 1, 0.5)

h = 0.0001
littlet = 1
bigt = 20

nsteps = ceiling(bigt/h)
nsaves = ceiling(bigt/littlet)
hilt = ceiling(littlet/h)
stopifnot((nsteps == (nsaves*hilt)))

h12 = sqrt(h)
xtraj = matrix(0, nrow = 1, ncol = (nsaves + 1))

# initial condition centered at origin so both -ve and +ve initial conditions generated
xtraj[,1] = 1

for (i in c(1:nsaves))
{
  # print loop counter
  print(i)
  flush.console()
  
  x = xtraj[,i]
  
  for (j in c(1:hilt))
    x = x + (thetavec[1])*(thetavec[2] - x)*(h) + (h12)*(thetavec[3])*(rnorm(n=1))      
  xtraj[,(i+1)] = x
}

tvec = seq(from = 0, to = bigt, by = littlet)
plot(tvec, xtraj, 'o')
xtraj = rbind(tvec, xtraj)
save(xtraj, file = 'fakedata.RData')
