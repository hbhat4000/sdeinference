driftfun <- function(c0, y, hermitefuncmat) 
{
  # computes drift as a vector
  matsize = dim(hermitefuncmat)
  mat = t(replicate(matsize[1],c0))*hermitefuncmat
  f = apply(mat,1,sum)
  return(f)
}

source('myhermfun.R')

#thetaapprox = c( -0.788892, -19.262578, -4.326613, -3.189801, -8.652582, 30.919872 ) #h=0.1
#thetaapprox = c( 0.989756, -30.419127, 7.786919, -70.145409, 6.537864, -36.557419 )#h=0.05
#thetaapprox = c( -1.177420, -61.319568, 1.818683, -134.410804, -0.181833, -31.680390 )#h=0.025
#thetaapprox = c( -1.552747, -78.447633, 1.576413, -172.979103, -0.599172, -40.674347 )#h=0.02
#thetaapprox = c( -0.363082, -9.178766, -1.669186, -20.471002, -1.568932, -9.106910 )#h=0.01
#thetaapprox = c ( -150.909020, -871.240082, -589.305878, -1298.032221, -967.594857, -505.911182 )#h=0.005
#thetaapprox = c( -5.616859, -261.057110, -1.792460, -550.225891 )#h=0.005

#thetaapprox = c( 0.010919, -5.417972, 0.163114, -11.367149 )#h=0.1
#thetaapprox = c( 0.070761, -4.643477, -0.090595, -9.396335 )#h=0.05
#thetaapprox = c( 0.074680, -4.171211, -0.126985, -8.424306 )#h=0.025
#thetaapprox = c( 0.073067, -3.873376, -0.136208, -7.811245 )#h=0.005
#thetaapprox =  c( 0.522027, -18.888363, 3.674679, -53.600928, 7.089392, -50.868739, 4.617129, -18.729303 )#h=0.005
#thetaapprox=c( 3.110242, -12.242620, 22.797122, -25.553821, 50.651185, -1.889225, 46.898711, 20.962576, 15.240726, 12.426160 )#h=0.005
#thetaapprox = c( -0.352049, -8.381179, -1.515531, -18.627244, -1.438252, -7.824235 )#H=0.005
#fakedata created with littlet=1,bigt=1, and 100 E-M samples



#thetaapprox=c( 0.114428, -28.165778, 0.196685, -64.866417, 0.236490, -35.211892 )#h=0.05
#thetaapprox=c( 0.019330, -10.239334, -0.206138, -22.793570, -0.136419, -10.504526 )#h=0.01
#thetaapprox=c( 0.025458, -9.664531, -0.180580, -21.452530, -0.112063, -9.681559 )  #h=0.005
thetaapprox=c( 0.261029, -12.968365, 0.926878, -26.695458, 0.724951, -2.495652, -0.990454, 20.129598, -1.143489, 11.359086 )#h=0.005


dx = 0.01
xvec = seq(-5,5,dx)
npts = length(xvec)
numpara = length(thetaapprox)
hermitefuncmat = myhermfun(numterms = numpara, xgrid = xvec)
fapprox = driftfun(thetaapprox, xvec, hermitefuncmat) 
ftrue = 4*(xvec - xvec^3)
plot(xvec, fapprox, type='l', col='red')
lines(xvec, ftrue )

error = sum((ftrue-fapprox)^2)*dx
print(error)







