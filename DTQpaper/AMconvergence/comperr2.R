# Solve the Example 3 SDE:
# dX_t = -(sin(X_t) cos^3(X_t)) dt + (cos^2(X_t)) dW_t
# with X_0 = 0

source('integrandmat.R')

comperr <- function(h,T,s,init,driftfun,difffun,exactfun)
{
    numsteps = ceiling(T/h)
    yM = pi/2 - (h^(2*s))
    bigM = ceiling(yM/(h^s))
    k = yM/bigM
    xvec = c(-bigM:bigM)*k
    npts = length(xvec)

    # pdf after one timestep
    mymu = init + driftfun(init)*h
    mysigma = abs(difffun(init))*sqrt(h)
    approxpdf = as.matrix(dnorm(x=xvec,mean=mymu,sd=mysigma))

    A = integrandmat(xvec,k,h,driftfun,difffun)
    for (i in c(1:npts)) A[,i] = A[,i] / (sum(A[,i])*k)

    # pdf after two timesteps
    for (i in c(2:numsteps))
        approxpdf = k*(A%*%approxpdf)

    truepdf = exactfun(xvec,T)

    pdferrinf = max(abs(approxpdf-truepdf))
    pdferrl1 = sum(abs(approxpdf-truepdf))*k
    ksnorm = max(abs(cumsum(approxpdf)*k-cumsum(truepdf)*k))

    errors = list(errinf=pdferrinf,errl1=pdferrl1,ks=ksnorm,yM=yM,npts=npts)
    return(errors)
}

# plot(xvec,approxpdf)
# lines(xvec,truepdf,col='red')


