rm(list = ls(all = TRUE))

# creating the runner's trajectory
T = 8

# setting the parameters for the pursuit model
speedchaser <- function(localt)
{
  if (localt < 4) sval = .4
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

h12 = sqrt(h)

xrun <- vector(length = nsteps + 1)
yrun <- vector(length = nsteps + 1)
xchase <- vector(length = nsteps + 1)
ychase <- vector(length = nsteps + 1)

for (i in c(0:nsteps))
{
  xrun[i+1] = -1 + 2*i/nsteps
  yrun[i+1] = sin(pi*xrun[i+1])
}
xchase[1] = -1 + rnorm(n=1,mean=0,sd=0.1)
ychase[1] = 1 + rnorm(n=1,mean=0, sd=0.1)

for (i in c(1:nsteps))
{
    # print loop counter
    xc = xchase[i]
    yc = ychase[i]
    xr = xrun[i]
    yr = yrun[i]
    
    localt = i*h
    xc = xc + (f1(xc, yc, xr, yr, localt))*(h) + (h12)*(g1(xc, yc, localt))*(rnorm(n = 1))
    yc = yc + (f2(xc, yc, xr, yr, localt))*(h) + (h12)*(g2(xc, yc, localt))*(rnorm(n = 1))

    xchase[i+1] = xc
    ychase[i+1] = yc
}

# make different versions of our one big simple path
ndiv = c(4000,2000,1000,500,400,200,100)

for (i in c(1:length(ndiv)))
{
    deltat = h*ndiv[i]
    myseq = seq(from=1,to=(nsteps+1),by=ndiv[i])
    lxrun = xrun[myseq]
    lyrun = yrun[myseq]
    lxchase = xchase[myseq]
    lychase = ychase[myseq]
    tvec = seq(from=0,to=T,by=deltat)
    runner = list(tvec, lxrun, lyrun)
    chaser = list(tvec, lxchase, lychase)
    myfname = paste('fakedataswitch4/fakedata_dt_',deltat,'.RData',sep='')
    save(runner, chaser, file=myfname)
}


