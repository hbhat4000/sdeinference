rm(list=ls(all=TRUE))

source('allcoeffs.R')

newhseq = c(0.005,0.002,0.001)
T = 1
s = 0.75
rootdir = "./cpp/ex"

for (j in c(1,2,4,5,6))
{
    newpdferrinf = numeric(length=3)
    newpdferrl1 = numeric(length=3)
    newksnorm = numeric(length=3)
    for (ell in c(1:3))
    {
        newfname = paste(rootdir,j,"/solution_",j,"_",newhseq[ell],sep='')
        numsteps = T/newhseq[ell]
        test = read.table(newfname,skip=(numsteps+2))
        approxpdf = test[,1]
        if (j==2) approxpdf = approxpdf[-length(approxpdf)]
        k = newhseq[ell]^s
        yM = pi/k
        bigm = ceiling(yM/k)
        xvec = seq(from=-bigm,to=bigm,by=1)*k
        truepdf = examples[[j]]$exact(xvec,T)
        print(c(length(approxpdf),length(truepdf)))
        newpdferrinf[ell] = max(abs(approxpdf-truepdf))
        newpdferrl1[ell] = sum(abs(approxpdf-truepdf))*k
        newksnorm[ell] = max(abs(cumsum(approxpdf)*k-cumsum(truepdf)*k))
    }
    fname = paste("newallresults_",j,".RData",sep="")
    save.image(fname)
}

