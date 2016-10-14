rm(list=ls(all=TRUE))

prefixes=c('./polydomain/','./logdomain/')
fnames=c('timingresults.RData','timingresultscpp.RData','timingresultssparse.RData')

alldata=list(NULL)
count=1
fplab=rep("FP",6)
lab=list(NULL)
lab[[1]]=rep("DTQ-Naive",6)
lab[[2]]=rep("DTQ-CPP",6)
lab[[3]]=rep("DTQ-Sparse",6)
domlab=list(NULL)
domlab[[1]]=rep("Linear",6)
domlab[[2]]=rep("Log",6)
for (iii in c(1:2))
{
  for (jjj in c(1:3))
  {
    thisfname=paste(prefixes[iii],fnames[jjj],sep='')
    load(thisfname)
    if (jjj==1)
    {
      alldata[[count]]=data.frame(error=FPl1errs/numreps,time=FPtimings/numreps,method=fplab,domain=domlab[[iii]])
      count = count + 1
    }
    alldata[[count]]=data.frame(error=DTQl1errs/numreps,time=DTQtimings/numreps,method=lab[[jjj]],domain=domlab[[iii]])
    count = count + 1
  }
}
finaldf = data.frame()
for (i in c(1:length(alldata)))
{
  finaldf = rbind(finaldf,alldata[[i]])
}
finaldf$method=as.factor(finaldf$method)
finaldf$domain=as.factor(finaldf$domain)

library('ggplot2')
library('grid')
library('gridExtra')

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl$nrow=1
  grid.arrange(do.call(arrangeGrob,gl),legend,nrow = 1,widths = unit.c(unit(1, "npc") - lwidth, lwidth))
}

myplots=list(NULL)

for (j in c(1:2))
{
  if (j==1)
    finaldflin=finaldf[which(finaldf$domain=="Linear"),]
  else
    finaldflin=finaldf[which(finaldf$domain=="Log"),]
  
  g = ggplot(finaldflin,aes(x=error,y=time,group=method))
  g = g + geom_point(aes(color=method)) + geom_line(aes(color=method))
  g = g + scale_x_continuous(trans='log10',breaks=c(0.003,0.01,0.03,0.10))
  if (j==1)
    g = g + scale_y_continuous(trans='log10',breaks=c(0.001,0.01,0.1,1,10,100))
  else
    g = g + scale_y_continuous(trans='log10',breaks=c(0.001,0.003,0.01,0.03,0.1,0.3))
  g = g + theme_bw()
  if (j==1)
    g = g + labs(title=expression(y[M] %prop% h^{-3/4}))
  else
    g = g + labs(title=expression(y[M] %prop% -log(h)))

  myplots[[j]] = g
}

grow = grid_arrange_shared_legend(myplots[[1]],myplots[[2]],ncol=2)
ggsave(filename="timings.pdf",plot=grow,height=4,width=8.5,units="in")
ggsave(filename="timings.eps",plot=grow,height=4,width=8.5,units="in")


