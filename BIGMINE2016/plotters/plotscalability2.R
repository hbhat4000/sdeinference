#clear memory
rm(list=ls(all=TRUE))
library('ggplot2')

#data
mylens=c(3,6,12,24)

#h=0.01
res1=numeric(length=4)
res1[1]=14+(4.519/60)
res1[2]=8+(11.491/60)
res1[3]=6+(1.702/60)
res1[4]=5+(15.841/60)

#h=0.02
res2=numeric(length=4)
res2[1]=3+(16.502/60)
res2[2]=2+(2.190/60)
res2[3]=1+(35.298/60)
res2[4]=1+(21.572/60)

alldata1=data.frame(L=mylens,time=res1,p="h=0.01")
alldata2=data.frame(L=mylens,time=res2,p="h=0.02")
alldata=rbind(alldata1,alldata2)

fit.1=lm(log(res1)~log(mylens))
fit.2=lm(log(res2)~log(mylens))
mylensart=seq(from=min(mylens),to=max(mylens),length.out=1000)
cur.1=exp(fit.1$coefficients[1])*mylensart^(fit.1$coefficients[2])
cur.2=exp(fit.2$coefficients[1])*mylensart^(fit.2$coefficients[2])
curdf.1=data.frame(L=mylensart,time=cur.1,p="h=0.01")
curdf.2=data.frame(L=mylensart,time=cur.2,p="h=0.02")

g = ggplot(alldata,aes(x=L,y=time,colour=p))
g = g + geom_point(size=2) 
g = g + xlab("number of Spark processors") + ylab("time (minutes)")
g = g + coord_trans(x="log10",y="log10")
g = g + geom_line(data=curdf.1)
g = g + geom_line(data=curdf.2)
g = g + scale_colour_manual("DTQ step",breaks=c("h=0.02","h=0.01"),values=c("blue","black"))
g = g + theme_bw()

ggsave(filename="scaling2.pdf",plot=g,height=3,width=4)
ggsave(filename="scaling2.eps",plot=g,height=3,width=4)


