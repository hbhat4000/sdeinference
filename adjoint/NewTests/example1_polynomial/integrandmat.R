integrandmat <- function(xvec, yvec, h, driftfun, difffun, c0)
{
  X = replicate(length(yvec), xvec)
  Y = t(replicate(length(xvec), yvec))
  thisdrift = driftfun(c0,Y)
  thisdiff = abs(difffun(c0, Y))
  out = exp(-(X - Y - thisdrift*h)^2/(2*thisdiff^2*h))/(thisdiff*sqrt(2*pi*h))
  return(out)
}