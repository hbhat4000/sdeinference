rm(list=ls(all=TRUE))

library(ggplot2)

dtqres = read.csv('./dtqou20.out',header=FALSE)
colnames(dtqres) = c("theta2","theta3","loglik")
v = ggplot(dtqres, aes(theta2, theta3, z = loglik))
v = v + geom_contour(binwidth = 100)
test = dtqres[which.max(dtqres[,3]),1:2]
v = v + geom_point(aes(x=test$theta2,y=test$theta3))
maxlllab = paste("(",format(test$theta2,digits=3),", ",format(test$theta3,digits=2),")",sep='')
v = v + annotate("text",x=test$theta2,y=test$theta3+0.02,label=maxlllab)
v = v + xlab(TeX("$\\theta_2$"))
v = v + ylab(TeX("$\\theta_3$"))
ggsave(filename="OUdtq20.pdf",plot=v,units="in",width=4,height=3)
ggsave(filename="OUdtq20.eps",plot=v,units="in",width=4,height=3)
print(v)

