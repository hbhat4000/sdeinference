rm(list=ls(all=TRUE))

library(ggplot2)
library(latex2exp)

Jvec = c(4, 8, 24)
Kvec = c(3, 4, 4)

for (i in c(1:3))
{
    fname = paste('llcf_',Jvec[i],'_',Kvec[i],'.RData',sep='')
    load(file=fname)

    asres = cbind(ppmat[,1:2],lltest)
    colnames(asres) = c("theta1","theta2","loglik")
    asres = data.frame(asres)

    v = ggplot(asres, aes(theta1, theta2, z = loglik))
    v = v + theme_bw()
    if (Jvec[i]==24)
        mybw = 1000
    else
        mybw = 100
    v = v + geom_contour(binwidth=mybw)
    v = v + xlab(TeX("$\\theta_1$"))
    v = v + ylab(TeX("$\\theta_2$"))
    basefname = paste('BIAS_',Jvec[i],'_',Kvec[i],sep='')
    pdffname = paste(basefname,'.pdf',sep='')
    epsfname = paste(basefname,'.eps',sep='')
    ggsave(filename=pdffname,plot=v,units="in",width=4,height=3)
    ggsave(filename=epsfname,plot=v,units="in",width=4,height=3)
}

