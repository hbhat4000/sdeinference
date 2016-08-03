rm(list = ls(all = TRUE))

for (iii in c(1:100))
{

# creating the runner's trajectory
t = 0.4
T = 0.4
tvec = seq(from = 0, to = T, by = t)

# xrun = seq(from = 1, to = xmax[1], by = 0.001)
# # yrun = 5*log(xrun) 
# yrun = abs(rnorm(n = length(xrun), mean = 0, sd = 2) + xrun*xrun/500 + xrun/20 + 1)
# xrun = xrun[seq(from = 1, to = length(xrun), length.out = T/t + 1)]
# yrun = yrun[seq(from = 1, to = length(yrun), length.out = T/t + 1)]

# setting the parameters for the pursuit model
speedchaser <- function(t)
{
  if (t <= 0.2) sval = .25
  else sval = .75
  return(sval)
}

f1 <- function(x, y, xr, yr, t)
{
  numerator = speedchaser(t) * (xr - x)
  denominator = sqrt((xr - x)^2 + (yr - y)^2)
  return(numerator/denominator)
}

f2 <- function(x, y, xr, yr, t)
{
  numerator = speedchaser(t) * (yr - y)
  denominator = sqrt((xr - x)^2 + (yr - y)^2)
  return(numerator/denominator)
}

g1 <- function(x, y, t)
{
  nu1 = 0.5
  return(nu1)
}

g2 <- function(x, y, t)
{
  nu2 = 0.5
  return(nu2)
}

# simulating the pursuit model
h = 0.4

nsteps = ceiling(T/h)   # 250000
nsaves = ceiling(T/t)   # 25
hilt = ceiling(t/h)		# 1000
stopifnot((nsteps == (nsaves*hilt)))

h12 = sqrt(h)

xrun <- vector(length = nsaves + 1)
yrun <- vector(length = nsaves + 1)
xchase <- vector(length = nsaves + 1)
ychase <- vector(length = nsaves + 1)

xrun[1] = abs(rnorm(n = 1, mean = 0, sd = 1)) 
yrun[1] = exp(xrun[1])
# yrun[1] = abs(rnorm(n = 1, mean = 0, sd = 1)) 
xchase[1] = abs(rnorm(n = 1, mean = xrun[1], sd = 0.5)) 
ychase[1] = abs(rnorm(n = 1, mean = yrun[1], sd = 0.5))

for (i in c(1:nsaves))
{
    # print loop counter
    print(i)
    flush.console()

    xc = xchase[i]
    yc = ychase[i]
    t = tvec[i]
    xr = xrun[i]
    yr = yrun[i]
    
    xrun[i+1] = abs(xrun[i] + rnorm(n = 1, mean = 0, sd = 0.25))
    yrun[i+1] = abs(yrun[i] + rnorm(n = 1, mean = 0, sd = 0.25))

    for (j in c(1:hilt))
    {
        xc = abs(xc + (f1(xc, yc, xr, yr, t))*(h) + (h12)*(g1(xc, yc, t))*(rnorm(n = 1)))     
        yc = abs(yc + (f2(xc, yc, xr, yr, t))*(h) + (h12)*(g2(xc, yc, t))*(rnorm(n = 1)))
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

myfname = paste('fakedata100/fakedata_h_0.4_',iii,'.RData',sep='')

save(runner, chaser, file=myfname)

}
