rm(list=ls(all=TRUE))

library('ggplot2')
library('scales')
library('grDevices')

plotlist = c(1,2,4,5,6)
for (i in c(1:5))
{
  j = plotlist[i]
  fname = paste("allresults_",j,".RData",sep="")
  load(fname)
  newfname = paste("newallresults_",j,".RData",sep="")
  load(newfname)

  fullhseq = hseq
  fullhseq = fullhseq[-length(fullhseq)]
  fullhseq = c(fullhseq, newhseq)

  fulll1errs = l1errs
  fulll1errs = fulll1errs[-length(fulll1errs)]
  fulll1errs = c(fulll1errs, newpdferrl1)

  linmod = lm(log(fulll1errs[-1]) ~ log(fullhseq[-1]))
  mypreds = data.frame(x=fullhseq[-1],y=as.numeric(exp(linmod$fitted.values)))

  fullinferrs = inferrs
  fullinferrs = fullinferrs[-length(fullinferrs)]
  fullinferrs = c(fullinferrs, newpdferrinf)

  fullkserrs = kserrs
  fullkserrs = fullkserrs[-length(fullkserrs)]
  fullkserrs = c(fullkserrs, newksnorm)
  
  myl1 = data.frame(x=fullhseq,y=fulll1errs)
  myinf = data.frame(x=fullhseq,y=fullinferrs)
  myks = data.frame(x=fullhseq,y=fullkserrs)

  mydat = data.frame(x=fullhseq,y=fulll1errs,norm="L1")
  mydat = rbind(mydat,data.frame(x=fullhseq,y=fullinferrs,norm="Linf"))
  mydat = rbind(mydat,data.frame(x=fullhseq,y=fullkserrs,norm="K-S"))

  myplot <- ggplot(data=mydat, aes(x=x,y=y,group=norm,color=norm))
  myplot <- myplot + theme_bw() + theme(plot.background = element_rect(fill='white'))

  myxticks = sort(10^(round(log(fullhseq)/log(10)*10)/10))
  rawyticks = round(log(mydat$y)/log(10)*10)/10
  rawyticks = round(seq(from=min(rawyticks),to=max(rawyticks),length.out=length(myxticks))*1)/1
  myyticks = unique(10^rawyticks)

  myplot <- myplot + scale_x_log10(breaks = fullhseq)
  myplot <- myplot + theme(axis.text.x = element_text(angle=90,hjust=1))

  myplot <- myplot + scale_y_log10(breaks = myyticks,
                                   labels = trans_format("log10", math_format(10^.x)))
  
  mytitle = paste("Example",j,sep=" ")
  # myplot <- myplot + labs(title = mytitle, x="h (temporal step size)", y="error")
  myplot <- myplot + labs(x="h (temporal step size)", y="error")
  myplot <- myplot + annotate("text",label=mytitle,y=max(l1errs),x=0.002)
  myplot <- myplot + geom_line() + geom_point()

  myplot <- myplot + geom_line(aes(x=x,y=y,group="L1"),data=mypreds,col="black")
  myslope = 0.001*round(1000*as.numeric(linmod$coefficients[2]))
  myplot <- myplot + annotate("text",label=paste("slope=",myslope,sep=''),x=fullhseq[6],y=fulll1errs[6]+0.01,size=4)

  fname = paste("convplot_",j,".eps",sep="")
  ggsave(filename=fname, plot=myplot, width=5, height=4, device=cairo_ps)

  fname = paste("convplot_",j,".pdf",sep="")
  ggsave(filename=fname, plot=myplot, width=5, height=4)
}


