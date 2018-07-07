source('integrandmat.R')

comperr <- function(h,T,s,init,driftfun,difffun,exactfun)
{
    numsteps = ceiling(T/h)
    k = h^s
    yM = k*(pi/(k^2))
    xvec = seq(-yM,yM,by=k)
    npts = length(xvec)

    # pdf after one timestep
    mymu = init + driftfun(init)*h
    mysigma = abs(difffun(init))*sqrt(h)
    approxpdf = as.matrix(dnorm(x=xvec,mean=mymu,sd=mysigma))

    A = integrandmat(xvec,xvec,h,driftfun,difffun)

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


