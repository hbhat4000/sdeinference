# Solve the SDE: dX_t = 1/2 X_t dt + sqrt(1 + X_t^2)dW_t
# with X_0 = 0

comperr <- function(h,T,s,init,driftfun,difffun,exactfun)
{
    ptm = proc.time()
    numsteps = ceiling(T/h)
    k = h^s
    yM = k*(pi/(k^2))
    xvec = seq(-yM,yM,by=k)
    npts = length(xvec)

    # pdf after one timestep
    fy <- function(y) exp(-(y-init-driftfun(init)*h)^2/(2*difffun(init)^2*h))/sqrt(2*pi*difffun(init)^2*h)

    A = k*integrandmat(xvec,xvec,h,driftfun,difffun)

    # pdf after two timesteps
    approxpdf = Matrix(fy(xvec),ncol=1)
    for (i in c(2:numsteps))
        approxpdf = A %*% approxpdf

    timetaken = proc.time() - ptm
    truepdf = exactfun(xvec,T)

    pdferrl1 = sum(abs(approxpdf-truepdf))*k
    return(list(timing=timetaken,error=pdferrl1))
}

