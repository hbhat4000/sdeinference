integrandmat <- function(xvec,yvec,h,driftfun,difffun)
{
  drifty = driftfun(yvec)
  diffy = difffun(yvec)
  mysd = diffy*sqrt(h)
  myvar = diffy*diffy*h
  
  c0 = 1/sqrt(2*pi)
  c0mod = c0/mysd
  dx = xvec[2] - xvec[1]
  bign = length(xvec)

  kvec = c(0)
  mydiags = list(NULL)

  # main diagonal
  maindiag = exp(-(h/2)*(drifty*drifty)/(diffy*diffy))*c0mod
  refsum = sum(abs(maindiag))*dx
  mydiags[[1]] = maindiag

  # super diagonals
  done = 0
  nk = 1
  mytol = 2e-16
  while (done==0)
  {
    curdiag = kvec[nk] + 1
    mymean = curdiag*dx + drifty*h
    thisdiag = exp(-mymean*mymean/(2*myvar))*c0mod
    thisdiag = thisdiag[-c(1:curdiag)]
    if ((nk==1) || (sum(abs(thisdiag)) > mytol*refsum))
    {
      mydiags[[nk+1]] = thisdiag
      kvec = c(kvec,curdiag)
      nk = nk + 1
    }
    else
    {
      done = 1
    }
  }
  kvec = c(kvec,-kvec[2:nk])
  for (i in c((nk+1):length(kvec)))
  {
    curdiag = kvec[i]
    mymean = curdiag*dx + drifty*h
    thisdiag = exp(-mymean*mymean/(2*myvar))*c0mod
    thisdiag = thisdiag[-c((bign+curdiag+1):bign)]
    mydiags[[i]] = thisdiag
  }
  
  out = bandSparse(bign,k=kvec,diagonals=mydiags)
  return(out)
}

