rm(list=ls(all=TRUE))

load('./run4/posteriorsamples_val4.RData')
x1 = x
den1 = density(x1^2)

load('./run5/posteriorsamples_val5.RData')
x2 = x
den2 = density(x2^2)

load('./run6/posteriorsamples_val6.RData')
x3 = x
den3 = density(x3^2)

xmin = min(c(den1$x,den2$x,den3$x))
xmax = max(c(den1$x,den2$x,den3$x))
ymin = min(c(den1$y,den2$y,den3$y))
ymax = max(c(den1$y,den2$y,den3$y))

postscript(file='densities.eps',width=5,height=5)

plot(den1$x,den1$y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type='l',col='blue',xlab='1/L',ylab='probability density')
lines(den2$x,den2$y,col='red')
lines(den3$x,den3$y,col='black')

empmode = den3$x[which.max(den3$y)]
abline(v=empmode)

abline(v=2*pi,lty=2)

legend(x="topright",legend=c("h=0.04","h=0.02","h=0.01"),col=c('blue','red','black'),lwd=1)

dev.off()
