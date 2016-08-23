myh = 0.001
myk = myh^0.75
mybigm = ceiling(pi/(myk^1.1))
load('fakedata.RData')
fd = xtraj
source('dtq_with_grad.R')
objgradfun <- function(c0)
{
    probmat = cdt(c0, h = myh, k = myk, bigm = mybigm, littlet = 1, data = fd)
    mylik = probmat$lik
    mylik[mylik < 0] = 0
    objective = -sum(log(mylik))
    nc0 = length(c0)
    gradient = numeric(length=nc0)
    for (i in c(1:nc0))
        gradient[i] = -sum(probmat$grad[[i]] / probmat$lik)

    return(list("objective"=objective,"gradient"=gradient))
}
objfun <- function(c0)
{
    return(objgradfun(c0)$objective)
}
gradfun <- function(c0)
{
    return(objgradfun(c0)$gradient)
}

optimHess(par=res$solution, fn=objfun, gr=gradfun)

