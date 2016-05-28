
driftfun <- function(c0, y) 
{
  # computes drift as a matrix
  nc0 = length(c0) - 1
  if (class(y)=="numeric") y = as.matrix(y) 
  f = matrix(data = c0[1],nrow=nrow(y),ncol=ncol(y))
  for (i in c(2:nc0))
    f = f + c0[i]*y^(i-1)

  # f = c0[1]*y + c0[2]*y^3
  return(f)
}

source('myhermfun.R')
#thetaapprox = c( 0.003736, 0.097670, -0.005797, -0.072356, 0.003621, 0.121049 ) #h=0.1
#thetaapprox = c( 0.440585, 1.032538, -0.184776, -0.869670, -0.125741, 0.124429 ) #h=0.1
#thetaapprox =c( 0.139511, 0.924743, 0.289763, 0.101952, 0.395997, -0.283494 )#h=0.05
#thetaapprox =  c( 0.015798, 2.106163, -0.353432, -0.174613, -0.029194, -1.079311, 0.422645, -0.293697, -0.122036, 0.641786 )#h=0.1
#thetaapprox = c( 0.200774, 6.540669, -0.413342, -2.841584, 0.260850, -3.628259, -0.249683, 0.381349, 0.052816, -0.580474 )#h=0.1
#above results are for the ata with littlet=1 and bigt=25,100 E-M samples


#thetaapprox = c( 0.071364, 2.517271, -0.037629, -2.457071, 0.683084 )#h=0.1
#thetaapprox = c( 0.054618, 3.395583, -0.010304, -3.356801, 0.837413 )#h=0.05
#thetaapprox = c( 0.062611, 2.518027, 0.002050, -2.457703, -0.026561, 0.683119 )#h=0.1
#thetaapprox =c ( 0.053826, 7.787202, 0.026081, -3.577966, -0.339615, -5.568593, 0.475490, 1.650997, -0.159996, 0.592318 )#h=0.1
#thetaapprox =c( 0.335946, 5.276142, 2.880895, 0.283317, -0.519049, -3.925081, -3.694777, -0.379138, 1.333054, 0.855038 )#h=0.05
#above results are for the ata with littlet=0.5 and bigt=25,



#thetaapprox = c( 0.029351, 2.347274, -0.024733, -2.293471, 0.667451 )#h=0.1
#thetaapprox = c( 3.755147, -0.679068, 4.257658, -2.561943, 188.790305 )#h=0.05
#thetaapprox =c( -1.384490, 3.922855, -1.589271, 3.419599, -3.354803, 4.337934, -5.617680, 171.353011 )#h=0.5
thetaapprox = c( 0.046674, 2.421273, -0.058239, -0.207580, -0.048950, -1.769787, 0.049930, 0.853305 )#h=0.05 and zero initial theta except diffcoeff=0.7
#above results are for the ata with littlet=0.5 and bigt=25,500 E-M samples


dx = 0.00001
xvec = seq(-5,5,dx)
npts = length(xvec)
numpara = length(thetaapprox)
fapprox = driftfun(thetaapprox, xvec) 
ftrue = 4*(xvec - xvec^3)
plot(xvec, fapprox, type='l', col='red')
lines(xvec, ftrue )