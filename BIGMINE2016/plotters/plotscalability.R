#clear memory
rm(list=ls(all=TRUE))
library('ggplot2')

#data
mylens=c(124,250,500,2500)

#h=0.01
res1=numeric(length=4)
res1[1]=(235+(30.585/60))/10
res1[2]=(439+(6.604/60))/10
res1[3]=88+(6.930/60)
res1[4]=419+(22.592/60)

#h=0.02
res2=numeric(length=4)
res2[1]=(59+(18.788/60))/10
res2[2]=(103+(25.361/60))/10
res2[3]=(197+(40.756/60))/10
res2[4]=89+(0.056/60)

alldata1=data.frame(L=mylens,time=res1,p="h=0.01")
alldata2=data.frame(L=mylens,time=res2,p="h=0.02")
alldata=rbind(alldata1,alldata2)

fit.1=lm(log(res1)~log(mylens))
fit.2=lm(log(res2)~log(mylens))
cur.1=exp(fit.1$coefficients[1])*mylens^(fit.1$coefficients[2])
cur.2=exp(fit.2$coefficients[1])*mylens^(fit.2$coefficients[2])
curdf.1=data.frame(L=mylens,time=cur.1,p="h=0.01")
curdf.2=data.frame(L=mylens,time=cur.2,p="h=0.02")

g = ggplot(alldata,aes(x=L,y=time,colour=p))
g = g + geom_point(size=2) 
g = g + xlab("L (length of observation series)") + ylab("time (minutes)")
g = g + coord_trans(x="log10",y="log10")
g = g + geom_line(data=curdf.1)
g = g + geom_line(data=curdf.2)
g = g + scale_colour_manual("DTQ step",breaks=c("h=0.02","h=0.01"),values=c("blue","black"))
g = g + theme_bw()
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename="scaling.pdf",plot=g,height=3,width=4)
ggsave(filename="scaling.eps",plot=g,height=3,width=4)

