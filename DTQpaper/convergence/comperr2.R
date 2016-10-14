# Solve the Example 3 SDE:
# dX_t = -(sin(X_t) cos^3(X_t)) dt + (cos^2(X_t)) dW_t
# with X_0 = 0

source('integrandmat.R')

comperr <- function(h,T,s,init,driftfun,difffun,exactfun)
{
    numsteps = ceiling(T/h)
    k = h^s
    yM = pi/2 - 2*k
    bigM = ceiling(yM/k)
    xvec = seq(from=-bigM,to=bigM,by=1)*k
    npts = length(xvec)

    # pdf after one timestep
    fy <- function(y) exp(-(y-init-driftfun(init)*h)^2/(2*difffun(init)^2*h))/sqrt(2*pi*difffun(init)^2*h)
    A = integrandmat(xvec,xvec,h,driftfun,difffun)
    for (i in c(1:npts)) A[,i] = A[,i] / (sum(A[,i])*k)

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


