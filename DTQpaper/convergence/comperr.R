# Solve the SDE: dX_t = 1/2 X_t dt + sqrt(1 + X_t^2)dW_t
# with X_0 = 0

source('integrandmat.R')

comperr <- function(h,T,s,init,driftfun,difffun,exactfun)
{
    numsteps = ceiling(T/h)
    k = h^s
    yM = k*(pi/(k^2))
    xvec = seq(-yM,yM,by=k)
    npts = length(xvec)

    # pdf after one timestep
    fy <- function(y) exp(-(y-init-driftfun(init)*h)^2/(2*difffun(init)^2*h))/sqrt(2*pi*difffun(init)^2*h)
    A = integrandmat(xvec,xvec,h,driftfun,difffun)

    # pdf after two timesteps
    approxpdf = k*(A%*%as.matrix(fy(xvec)))
    for (i in c(3:numsteps))
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


