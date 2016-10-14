rm(list=ls(all=TRUE))

library('Matrix')

source('integrandmatcpp.R')
source('comperr.R')
source('fpwithg.R')
source('allcoeffs.R')

hseq = c(0.5,0.2,0.1,0.05,0.02,0.01)
init = c(0,0,0,0,0,0,0,0,0)
numh = length(hseq)
T = 1
s = 0.75
# enums = c(3,2,8,5,4)
enums = c(3)
DTQl1errs = matrix(0,nrow=numh,ncol=length(enums))
DTQtimings = matrix(0,nrow=numh,ncol=length(enums))
# FPl1errs = matrix(0,nrow=numh,ncol=length(enums))
# FPtimings = matrix(0,nrow=numh,ncol=length(enums))

numreps = 100

for (rep in c(1:numreps))
{
  print(paste("Rep",rep))
  flush.console()
  for (jj in c(1:length(enums)))
  {
    j = enums[jj]
    print(paste("Example",j))
    flush.console()
    for (ii in c(1:numh))
    {
      DTQerr = comperr(hseq[ii],T,s,init[j],drift=examples[[j]]$drift,diff=examples[[j]]$diff,exact=examples[[j]]$exact)
      DTQl1errs[ii,jj] = DTQl1errs[ii,jj] + DTQerr$error
      DTQtimings[ii,jj] = DTQtimings[ii,jj] + as.numeric(DTQerr$timing[3])
      
      # FPerr = fperr(hseq[ii],T,s,drift=examples[[j]]$drift,diff=examples[[j]]$diff,exact=examples[[j]]$exact)
      # FPl1errs[ii,jj] = FPl1errs[ii,jj] + FPerr$error
      # FPtimings[ii,jj] = FPtimings[ii,jj] + as.numeric(FPerr$timing[3])
    }
    save.image(file="timingresultscpp.RData")
  }
}

