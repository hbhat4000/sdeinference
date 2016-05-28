#clear memory
rm(list=ls(all=TRUE))
library('ggplot2')
library('reshape2')

#pick which directories and which variables we want
mydirs=c(7)
myvbls=c(4:254)

#colors
mycols=c("black","blue")
mypars=c("h=0.02","h=0.01")

#set burn-in period
burnin=c(1:100)

#load data
basedir="actual_bigrun"
j=1
jj=mydirs[j]
fname=paste(basedir,jj,"/mcmc.out",sep='')
temp=read.csv(file=fname,header=FALSE)
mydata=temp[,myvbls]
sigeps=sqrt(mean(temp[,255]))
mydata=mydata[-burnin,]

#add time vector
fname=paste(basedir,jj,"/tvec.csv",sep='')
tvec=read.csv(file=fname,header=FALSE)
names(tvec)=names(mydata)
mydata=rbind(tvec,mydata)

#load yvec vector (yes, it's called xvec)
fname=paste(basedir,jj,"/xvec.csv",sep='')
yvec=read.csv(file=fname,header=FALSE)
ydf=data.frame(tval=as.numeric(tvec),value=as.numeric(yvec))
ydf$variable="y"

#transpose and melt
mydata=data.frame(t(mydata))
colnames(mydata)[1]="tval"
mydf = melt(mydata,id.vars="tval")

#calculate mean(xvec)
xvec=as.numeric(rowMeans(mydata[,-1]))
xdf=data.frame(tval=as.numeric(tvec),value=xvec)
xdf$variable="x"

g = ggplot(mydf,aes(x=tval,y=value,group=variable)) + geom_line(size=0.2,colour="grey50") 
g = g + geom_line(data=ydf,colour="red",size=0.5)
g = g + geom_line(data=xdf,colour="black",size=0.5)
g = g + xlab("t") + ylab("observations and inferred states")
g = g + theme_bw()

ggsave(filename="timeseries1.pdf",plot=g,height=3,width=6,units="in")
ggsave(filename="timeseries1.eps",plot=g,height=3,width=6,units="in")

ydf$lower = ydf$value-sigeps
ydf$upper = ydf$value+sigeps

g2 = ggplot(ydf,aes(x=tval,y=value)) + geom_line(colour="red",size=0.5)
g2 = g2 + geom_errorbar(aes(ymin=lower,ymax=upper),size=0.2,width=0.2,colour="grey30")
g2 = g2 + geom_line(data=xdf,colour="black",size=0.5)
g2 = g2 + theme_bw()
g2 = g2 + xlab("t") + ylab("observations with error bars,\n and inferred state")

ggsave(filename="timeseries2.pdf",plot=g2,height=3,width=6,units="in")
ggsave(filename="timeseries2.eps",plot=g2,height=3,width=6,units="in")







