# function that computes log likelihood
mylikelihood <- function(likden, dat)
{
    numsteps = dim(likden$y)[2]
    loglik = 0
    for (k in c(1:numsteps))
    {
        likdat = approx(x = likden$x, y = likden$y[,k], xout = dat[2:nrow(dat),k+1], yleft = 2.2e-16, yright = 2.2e-16)
        loglik = sum(log(likdat$y)) + loglik	
    }
    return(loglik)
}

gk <- function(likden, dat)
{
    numsteps = dim(likden$z1)[2]
    val1 = 0
    val2 = 0
    for (k in c(1:numsteps))
    {
       # using interpolation on P and gradients of P
       prob = approx(x = likden$x, y = likden$y[,k],  xout = dat[2:nrow(dat),k+1], yleft = 2.2e-16, yright = 2.2e-16)
       diffprob1 = approx(x = likden$x, y = likden$z1[,k], xout = dat[2:nrow(dat),k+1], yleft = 2.2e-16, yright = 2.2e-16)
       diffprob2 = approx(x = likden$x, y = likden$z2[,k], xout = dat[2:nrow(dat),k+1], yleft = 2.2e-16, yright = 2.2e-16)
       val1 = sum(diffprob1$y / prob$y) + val1
       val2 = sum(diffprob2$y / prob$y) + val2
    }
    
    vec = c(val1, val2)
    return(vec)
}

