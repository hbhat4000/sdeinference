rm(list = ls(all = TRUE))

for (iii in c(1:100))
{

# creating the runner's trajectory
deltat = 0.1
T = 8
tvec = seq(from = 0, to = T, by = deltat)

# xrun = seq(from = 1, to = xmax[1], by = 0.001)
# # yrun = 5*log(xrun) 
# yrun = abs(rnorm(n = length(xrun), mean = 0, sd = 2) + xrun*xrun/500 + xrun/20 + 1)
# xrun = xrun[seq(from = 1, to = length(xrun), length.out = T/t + 1)]
# yrun = yrun[seq(from = 1, to = length(yrun), length.out = T/t + 1)]

# setting the parameters for the pursuit model
speedchaser <- function(localt)
{
  if (localt < 4) sval = .2
  else sval = 1.0
  return(sval)
}

f1 <- function(x, y, xr, yr, localt)
{
  numerator = speedchaser(localt) * (xr - x)
  denominator = sqrt((xr - x)^2 + (yr - y)^2)
  return(numerator/denominator)
}

f2 <- function(x, y, xr, yr, localt)
{
  numerator = speedchaser(localt) * (yr - y)
  denominator = sqrt((xr - x)^2 + (yr - y)^2)
  return(numerator/denominator)
}

g1 <- function(x, y, localt)
{
  nu1 = 0.15
  return(nu1)
}

g2 <- function(x, y, localt)
{
  nu2 = 0.1
  return(nu2)
}

# simulating the pursuit model
h = 0.0001

nsteps = ceiling(T/h)        # 250000
nsaves = ceiling(T/deltat)   # 25
hilt = ceiling(deltat/h)     # 1000
stopifnot((nsteps == (nsaves*hilt)))

h12 = sqrt(h)

xrun <- vector(length = nsaves + 1)
yrun <- vector(length = nsaves + 1)
xchase <- vector(length = nsaves + 1)
ychase <- vector(length = nsaves + 1)

# xrun[1] = abs(rnorm(n = 1, mean = 0, sd = 1)) 
# yrun[1] = exp(xrun[1])
# yrun[1] = abs(rnorm(n = 1, mean = 0, sd = 1)) 
# xchase[1] = abs(rnorm(n = 1, mean = xrun[1], sd = 0.5)) 
# ychase[1] = abs(rnorm(n = 1, mean = yrun[1], sd = 0.5))
for (i in c(0:nsaves))
{
  xrun[i+1] = -1 + 2*i/nsaves
  yrun[i+1] = sin(pi*xrun[i+1])
}
xchase[1] = -1 + rnorm(n=1,mean=0,sd=0.1)
ychase[1] = 1 + rnorm(n=1,mean=0, sd=0.1)

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    xc = xchase[i]
    yc = ychase[i]
    xr = xrun[i]
    yr = yrun[i]
    
#    xrun[i+1] = abs(xrun[i] + rnorm(n = 1, mean = 0, sd = 0.25))
#    yrun[i+1] = abs(yrun[i] + rnorm(n = 1, mean = 0, sd = 0.25))

    for (j in c(1:hilt))
    {
        localt = tvec[i] + (j-1)*h
        xc = xc + (f1(xc, yc, xr, yr, localt))*(h) + (h12)*(g1(xc, yc, localt))*(rnorm(n = 1))
        yc = yc + (f2(xc, yc, xr, yr, localt))*(h) + (h12)*(g2(xc, yc, localt))*(rnorm(n = 1))
    }

    xchase[i+1] = xc
    ychase[i+1] = yc
}

xaxis = c(min(xchase, xrun), max(xchase, xrun))
yaxis = c(min(ychase, yrun), max(ychase, yrun))

runner = list(tvec, xrun, yrun)
chaser = list(tvec, xchase, ychase)

# plot(xrun, yrun, xaxis, yaxis, type = "l", col = "red")
# lines(xchase, ychase, xaxis, yaxis, type = "b", col = "black")

myfname = paste('fakedataswitch3/fakedata_h_0.0001_',iii,'.RData',sep='')

save(runner, chaser, file=myfname)

}
