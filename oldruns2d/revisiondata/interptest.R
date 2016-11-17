library('Rdtq2d')
options(digits=16)
X = matrix(0, nrow = 2, ncol = 3)
myh = 0.05
myns = 4
myk = 0.1
xylimit = 0.5
M = ceiling(xylimit/myk)
xvec = myk*c(-M:M)

X[1,] = c(0.1,0.2,0.0)
X[2,] = c(0.3,0.4,0.1)
mm = length(xvec)

source('truethetavec.R')
thetavec = truethetavec

interpolate = 1

# for "with interpolation"
if (interpolate==1) {
library('fields')
oldden = Rdtq2d(thetavec, as.matrix(X[1,1]), as.matrix(X[1,2]), myh, myns, myk, xylimit)
myposterior <- function(likden, dat)
{
    steps = ncol(likden)
    # for scenario 2, likden has 'datapoints' pdfs  
    logpost = 0
    for(i in c(1:steps))
    { 
      myloc = matrix(dat[(i+1),c(1:2)],ncol=2)
      likdat = interp.surface(obj=list(x=xvec,y=xvec,z=matrix(likden[,i],nrow=mm)),loc=myloc)
      likdat[likdat <= 2.2e-16] = 2.2e-16
      logpost = logpost + sum(log(likdat))
    }
    return(logpost)
}
finalden = myposterior(oldden, X)
}
if (interpolate==0) {
# for "no interpolation"
oldden = Rdtq2d(thetavec, X[,1], X[,2], myh, myns, myk, xylimit)
myposterior <- function(den)
{
  den[den <= 2.2e-16] = 2.2e-16
  loglik = sum(log(den))
  return(loglik)
}
finalden = myposterior(oldden)
}


##########
##### with 2 steps
# -4.388115943043199 : with interp, with gamma
# -4.387483267742011 : with interp, with gamma = veclen
# -4.392433102799039 : with interp, no gamma
# -4.392433102799038 : no interp

#### with 3 steps
# -5.032223868172064 : with interp, no gamma
# -5.032465647337961 : no interp

#### with 4 steps
#
# -5.838741048012646 : no interp
