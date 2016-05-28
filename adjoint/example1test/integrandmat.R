integrandmat <- function(xvec, yvec, h, driftfun, difffun, c0, hermitefuncmat, diffcoeff)
{
  X = replicate(length(yvec), xvec)
  Y = t(replicate(length(xvec), yvec))
  thisdriftvec = driftfun(c0,yvec,hermitefuncmat)
  thisdrift = t(replicate(length(xvec), thisdriftvec))
  thisdiff = abs(difffun(diffcoeff, Y))
  out = exp(-(X - Y - thisdrift*h)^2/(2*thisdiff^2*h))/(thisdiff*sqrt(2*pi*h))
  return(out)
}