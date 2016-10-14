integrandmat <- function(xvec,yvec,h,driftfun,difffun)
{
  X=replicate(length(yvec),xvec)
  Y=t(replicate(length(xvec),yvec))
  out = exp(-(X-Y-driftfun(Y)*h)^2/(2*difffun(Y)^2*h))/(difffun(Y)*sqrt(2*pi*h))
  return(out)
}
