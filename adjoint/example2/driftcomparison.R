driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

source('myhermfun.R')

# thetaapprox = c( -0.059761, 2.070630, 0.030253, 0.643363, -0.037453 )
# thetaapprox =  c( -0.086554, 2.058772, 0.024500, 1.125249, -0.043199, -0.640039 )
# thetaapprox = c( -5.829304, 28.725348, -36.354984, 12.031942, -18.446729, 9.577597, 39.282653, -56.016274 )
# Above three results are using myh=0.1 and initial theta = c(1:n)


# thetaapprox = c( -0.033477, 1.828456, 0.048299, 0.754616, -0.093856, 0.135999, 0.052135, -0.996307 ) # h=0.1 and random 8 init parameters
thetaapprox = c( -0.059957, 1.934193, 0.061708, 0.813082, -0.079857, -0.360314, 0.010801, -0.445946, 0.046713, -0.857468, -0.064780, -0.725007 )


dx = 0.0001
xvec = seq(-5,5,dx)
npts = length(xvec)
numpara = length(thetaapprox)
hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)
fapprox = driftfun(thetaapprox, xvec, hermitefuncmat) 
ftrue = sin(xvec)
plot(xvec, fapprox, type='l', col='red')
lines(xvec, ftrue )

error = sum((ftrue-fapprox)^2)*dx
print(error)

#[1] 2.442744 error with 5 hermite functions
#[1] 2.875677 error with 6 hermite functions
#[1] 7264.721 error with 6 hermite functions
# Above three results are using myh=0.1 and initial theta = c(1:n)


#[1] 2.054157 error with random 8 parameters for h=0.1
#[1] 0.658039 error with random 12 parameters for h=0.1