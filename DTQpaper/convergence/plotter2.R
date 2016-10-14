rm(list=ls(all=TRUE))

library('ggplot2')
library('scales')
library('grDevices')

for (j in c(3:3))
{
  fname = paste("allresults_",j,".RData",sep="")
  load(fname)

  fullhseq = hseq
  fulll1errs = l1errs

  myskip = c(1)
  linmod = lm(log(fulll1errs[-myskip]) ~ log(fullhseq[-myskip]))
  mypreds = data.frame(x=fullhseq[-myskip],y=as.numeric(exp(linmod$fitted.values)))

  fullinferrs = inferrs
  fullkserrs = kserrs
  
  myl1 = data.frame(x=fullhseq,y=fulll1errs)
  myinf = data.frame(x=fullhseq,y=fullinferrs)
  myks = data.frame(x=fullhseq,y=fullkserrs)

  mydat = data.frame(x=fullhseq,y=fulll1errs,norm="L1")
  mydat = rbind(mydat,data.frame(x=fullhseq,y=fullinferrs,norm="Linf"))
  mydat = rbind(mydat,data.frame(x=fullhseq,y=fullkserrs,norm="K-S"))

  myplot <- ggplot(data=mydat, aes(x=x,y=y,group=norm,color=norm))
  myplot <- myplot + theme_bw()

  myxticks = sort(10^(round(log(fullhseq)/log(10)*10)/10))
  rawyticks = round(log(mydat$y)/log(10)*10)/10
  rawyticks = round(seq(from=min(rawyticks),to=max(rawyticks),length.out=length(myxticks))*1)/1
  myyticks = unique(10^rawyticks)

  myplot <- myplot + scale_x_log10(breaks = fullhseq)
  myplot <- myplot + theme(axis.text.x = element_text(angle=90,hjust=1))

  myplot <- myplot + scale_y_log10(breaks = myyticks,
                                   labels = trans_format("log10", math_format(10^.x)))
  
  k = 3
  mytitle = paste("Example",k,sep=" ")
  # myplot <- myplot + labs(title = mytitle, x="h (temporal step size)", y="error")
  myplot <- myplot + labs(x="h (temporal step size)", y="error")
  myplot <- myplot + annotate("text",label=mytitle,y=max(l1errs),x=0.002)
  myplot <- myplot + geom_line() + geom_point()

  myplot <- myplot + geom_line(aes(x=x,y=y,group="L1"),data=mypreds,col="black")
  myslope = 0.001*round(1000*as.numeric(linmod$coefficients[2]))
  myplot <- myplot + annotate("text",label=paste("slope=",myslope,sep=''),x=fullhseq[5],y=fulll1errs[5]+0.01,size=4)

  fname = paste("convplot_",k,".eps",sep="")
  ggsave(filename=fname, plot=myplot, width=5, height=4, device=cairo_ps)

  fname = paste("convplot_",k,".pdf",sep="")
  ggsave(filename=fname, plot=myplot, width=5, height=4)
}


