rm(list=ls(all=TRUE))

library(ggplot2)
library(latex2exp)

dtqres = read.csv('./dtqbistab.out',header=FALSE)
colnames(dtqres) = c("theta1","theta2","loglik")
v = ggplot(dtqres, aes(theta1, theta2, z = loglik))
v = v + theme_bw()
v = v + geom_contour(binwidth = 100)
test = dtqres[which.max(dtqres[,3]),1:2]
v = v + geom_point(aes(x=test$theta1,y=test$theta2))
maxlllab = paste("(",format(test$theta1,digits=3),", ",format(test$theta2,digits=2),")",sep='')
v = v + annotate("text",x=test$theta1,y=test$theta2+0.2,label=maxlllab)
v = v + xlab(TeX("$\\theta_1$"))
v = v + ylab(TeX("$\\theta_2$"))
ggsave(filename="./BIdtqzoom.pdf",plot=v,units="in",width=4,height=3)
ggsave(filename="./BIdtqzoom.eps",plot=v,units="in",width=4,height=3)
print(v)


