rm(list=ls(all=TRUE))

source('comperr.R')
source('allcoeffs.R')

hseq = c(0.5,0.2,0.1,0.05,0.02,0.01)
init = c(0,0,0,0,0,0,0,0,0)
numh = length(hseq)
T = 1
s = 0.75
inferrs = numeric(length=numh)
l1errs = numeric(length=numh)
kserrs = numeric(length=numh)
yM = numeric(length=numh)
npts = numeric(length=numh)

#for (j in c(1,2,4,5,6))
for (j in c(5,6))
{
  for (ii in c(1:numh))
  {
    print(paste0("j=",j,", h=",hseq[ii]))
    err = comperr(hseq[ii],T,s,init[j],drift=examples[[j]]$drift,diff=examples[[j]]$diff,exact=examples[[j]]$exact)
    inferrs[ii] = err$errinf
    l1errs[ii] = err$errl1
    kserrs[ii] = err$ks
    yM[ii] = err$yM
    npts[ii] = err$npts
  }
  fname = paste("allresults_",j,".RData",sep="")
  save.image(file=fname)
  fname = paste("l1plot_",j,".eps",sep="")
  postscript(file=fname,width=5,height=5)
  plot(log(hseq),log(l1errs),type='o',col='blue',ann=FALSE)
  abline(lm(log(l1errs) ~ log(hseq)))
  title(xlab="log(h)")
  dev.off()
}


