integrandmat <- function(xvec,yvec,h,driftfun,difffun,c0)
{
  X=replicate(length(yvec),xvec)
  Y=t(replicate(length(xvec),yvec))
  thisdiff = abs(difffun(c0,Y))
  out = exp(-(X-Y-driftfun(c0,Y)*h)^2/(2*thisdiff^2*h))/(thisdiff*sqrt(2*pi*h))
  return(out)
}
