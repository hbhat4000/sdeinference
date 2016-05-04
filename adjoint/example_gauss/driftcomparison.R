driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

source('myhermfun.R')
#thetaapprox = c( 8.206340, -6.041504, 3.685637, -2.528507, -7.456732, 3.631032 )#h=0.1
#thetaapprox = c( 25.944129, -9.710085, -0.333961, -14.580044, 14.669668, -48.350622, 19.880589, -33.713189, -19.377095, 4.652954 )#h=0.1
# fakedata created with littlet=1, bigt=25, 100 E-M samples




#thetaapprox = c( 9.485894, -9.986393, 4.861590, 1.518902, -5.638759, -4.389839 )#h=0.1
#thetaapprox = c( 7.605415, -7.480514, 4.323765, 1.114767, -5.848125, 4.154041 )#h=0.05
#thetaapprox =  c( 8.428031, -9.644130, 2.519714, 6.087432, -4.144822, -4.829452, 0.141409, 4.161184 )#h=0.5
#thetaapprox = c( 5.553570, -4.250274, 0.321183, 1.858188, 0.674846, -3.188573, -2.857642, 3.591634 )#h=0.025


thetaapprox = c( 0.571901, 0.054604, 0.033615, 0.159754, -0.006645, -0.004924, -0.018042, 0.981503 )#h=0.005
#thetaapprox = c ( 0.577116, 0.055502, 0.039918, 0.162266, 0.980678 )#h=0.005
# fakedata created with littlet=0.5, bigt=25, 100 E-M samples
# zero init theta for drift coef and diff coef is 0.7



dx = 0.001
xvec = seq(-5,5,dx)
npts = length(xvec)
numpara = length(thetaapprox)-1
hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)
fapprox = driftfun(c0=thetaapprox[1:numpara], xvec, hermitefuncmat)
mu = 0
sigma = 1
ftrue = exp(-(xvec-mu)^2/(2*sigma^2))/sqrt(2*pi*sigma^2)
plot(xvec, fapprox, type='l', col='red')
lines(xvec, ftrue )
error = sum((ftrue-fapprox)^2)*dx
print(error)