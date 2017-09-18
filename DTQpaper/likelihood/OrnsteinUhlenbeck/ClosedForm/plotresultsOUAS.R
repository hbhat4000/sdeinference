rm(list=ls(all=TRUE))

library(ggplot2)
library(latex2exp)

load('./OUAS4_3.RData')
n0 = 50
n1 = 50
theta0grid = seq(from=1,to=7,length.out=n0)
theta1grid = seq(from=0.125,to=0.5,length.out=n1)
asres = matrix(0,nrow=n0*n1,ncol=3)
cr = 1
for (i0 in c(1:n0))
{
  for (i1 in c(1:n1))
  {
    asres[cr,1] = theta0grid[i0]
    asres[cr,2] = theta1grid[i1]
    asres[cr,3] = llmat[i0,i1]
    cr = cr + 1
  }
}
colnames(asres) = c("theta2","theta3","loglik")
asres = data.frame(asres)

v = ggplot(asres, aes(theta2, theta3, z = loglik))
v = v + theme_bw()
v = v + geom_contour(binwidth=100)
test = asres[which.max(asres[,3]),]
v = v + geom_point(aes(x=test$theta2,y=test$theta3))
#v = v + annotate("text",x=1.16,y=8.35,label="loglik = -443")
maxlllab = paste("(",format(test$theta2,digits=3),", ",format(test$theta3,digits=2),")",sep='')
v = v + annotate("text",x=test$theta2,y=test$theta3+0.02,label=maxlllab)
v = v + xlab(TeX("$\\theta_2$"))
v = v + ylab(TeX("$\\theta_3$"))
ggsave(filename="OUAS.pdf",plot=v,units="in",width=4,height=3)
ggsave(filename="OUAS.eps",plot=v,units="in",width=4,height=3)

