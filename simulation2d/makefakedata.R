rm(list = ls(all = TRUE))

# creating the runner's trajectory
t = 1
T = 100
tvec = seq(from = 0, to = T, by = t)
xrun = seq(from = 1, to = 94, by = 0.001)
# yrun = 5*log(xrun) 
yrun = rnorm(n = length(xrun), mean = 0, sd = 2) + xrun*xrun/500 + xrun/20 + 5
xrun = xrun[seq(from = 1, to = length(xrun), length.out = T+t)]
yrun = yrun[seq(from = 1, to = length(yrun), length.out = T+t)]

# Basketbal court's dimension
xaxis = c(0,94)
yaxis = c(0,50)
# plot(xrun, yrun, xaxis, yaxis, type = "o", col = "red")
plot(xrun, yrun, xaxis, yaxis, type = "l", col = "red")
# lines(xrun, yrun, type = "l", col = "blue")
runner = list(tvec, xrun, yrun)

# setting the parameters for the pursuit model
speed <- function(t)
{
  # sval = 1 + 2*t
  # sval = -1 + 2*t + rnorm(n = length(t))
  sval = 1.2 + rnorm(n = length(t))
  return(sval)
}

f1 <- function(x, y, xr, yr, t)
{
  numerator = speed(t) * (xr - x)
  denominator = sqrt((xr - x)^2 + (yr - y)^2)
  return(numerator/denominator)
}

f2 <- function(x, y, xr, yr, t)
{
  numerator = speed(t) * (yr - y)
  denominator = sqrt((xr - x)^2 + (yr - y)^2)
  return(numerator/denominator)
}

g1 <- function(x, y, t)
{
  nu1 = rnorm(n = 1, mean = 0.5, sd = 1)
  return(nu1)
}

g2 <- function(x, y, t)
{
  nu2 = rnorm(n = 1, mean = 0.5, sd = 1)
  return(nu2)
}

# simulating the pursuit model
h = 0.0001


nsteps = ceiling(T/h)   # 250000
nsaves = ceiling(T/t)   # 25
hilt = ceiling(t/h)		# 1000
stopifnot((nsteps == (nsaves*hilt)))

h12 = sqrt(h)

# xchase = matrix(0, nrow = ntrials, ncol = (nsaves + 1))
# xchase[,1] = rnorm(n = ntrials, mean = xrun[1], sd = sqrt(h)) 
# ychase = matrix(0, nrow = ntrials, ncol = (nsaves + 1))
# ychase[,1] = rnorm(n = ntrials, mean = yrun[1], sd = sqrt(h)) 

xchase <- vector(length = nsaves + 1)
ychase <- vector(length = nsaves + 1)
xchase[1] = rnorm(n = 1, mean = xrun[1], sd = sqrt(h)) 
ychase[1] = rnorm(n = 1, mean = yrun[1], sd = sqrt(h))

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
    
    for (j in c(1:hilt))
    {
        xc = xc + (f1(xc, yc, xr, yr, t))*(h) + (h12)*(g1(xc, yc, t))*(rnorm(n = 1))      
        yc = yc + (f2(xc, yc, xr, yr, t))*(h) + (h12)*(g2(xc, yc, t))*(rnorm(n = 1))
    }

    xchase[i + 1] = xc
    ychase[i + 1] = yc
}

lines(xchase, ychase, type = "b", col = "black")
chaser = list(tvec, xchase, ychase)
save(runner, chaser, file = 'fakedata.RData')