# rm(list=ls(all=TRUE))

# load('samples_pomp_vanderpol_by2.RData')
# x2 = mcmcsamples[,2]
# den2 = density(x2^2, bw = 'SJ')

# load('samples_pomp_vanderpol_by4.RData')
# x4 = mcmcsamples[,2]
# den4 = density(x4^2, bw = 'SJ')

# load('samples_pomp_vanderpol_by8.RData')
# x8 = mcmcsamples[,2]
# den8 = density(x8^2, bw = 'SJ')

# # xmin = min(c(den2$x,den4$x,den8$x))
# # xmax = max(c(den2$x,den4$x,den8$x))
# ymin = min(c(den2$y,den4$y,den8$y))
# ymax = max(c(den2$y,den4$y,den8$y))
# xmin = 0.5
# xmax = 1.5

# postscript(file='pomp_theta2_SJ.eps',width=5,height=5)
# # pdf(file='pomp_theta2_SJ.pdf',width=5,height=5)

# plot(den2$x,den2$y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type='l',col='blue',xlab='theta_2',ylab='probability density')
# lines(den4$x,den4$y,col='red')
# lines(den8$x,den8$y,col='black')

# empmode = den8$x[which.max(den8$y)]
# abline(v=empmode)

# abline(v=1,lty=2)

# legend(x="topright",legend=c("h=0.05","h=0.025","h=0.0125"),col=c('blue','red','black'),lwd=1)

# dev.off()

rm(list=ls(all=TRUE))

load('samples_vanderpol_by2.RData')
x2 = mcmcsamples[,3]
den2 = density(x2^2, bw = 'SJ')

load('samples_vanderpol_by4.RData')
x4 = mcmcsamples[,3]
den4 = density(x4^2, bw = 'SJ')

load('samples_vanderpol_by8.RData')
x8 = mcmcsamples[,3]
den8 = density(x8^2, bw = 'SJ')

# xmin = min(c(den2$x,den4$x,den8$x))
# xmax = max(c(den2$x,den4$x,den8$x))
ymin = min(c(den2$y,den4$y,den8$y))
ymax = max(c(den2$y,den4$y,den8$y))
xmin = 3.9
xmax = 4.7

# postscript(file='dtq_theta1.eps',width=5,height=5)
pdf(file='dtq_theta3_SJ.pdf',width=5,height=5)

plot(den2$x,den2$y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type='l',col='blue',xlab='theta_3',ylab='probability density')
lines(den4$x,den4$y,col='red')
lines(den8$x,den8$y,col='black')

empmode = den8$x[which.max(den8$y)]
abline(v=empmode)

abline(v=4,lty=2)

legend(x="topright",legend=c("h=0.05","h=0.025","h=0.0125"),col=c('blue','red','black'),lwd=1)

dev.off()
