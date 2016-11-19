rm(list=ls(all=TRUE))

load('./samples_by100/samples_by100_hby1.RData')
x1 = x
den1 = density(x1^2)

load('./samples_by100/samples_by100_hby2.RData')
x2 = x
den2 = density(x2^2)

# load('./samples_by100/samples_by100_hby3.RData')
# x3 = x
# den3 = density(x3^2)

# load('./samples_by100/samples_by100_hby4.RData')
# x4 = x
# den4 = density(x4^2)

# load('./samples_by100/samples_by100_hby5.RData')
# x5 = x
# den5 = density(x5^2)


# xmin = min(c(den1$x,den2$x,den3$x,den4$x,den5$x))
# xmax = max(c(den1$x,den2$x,den3$x,den4$x,den5$x))
# ymin = min(c(den1$y,den2$y,den3$y,den4$y,den5$y))
# ymax = max(c(den1$y,den2$y,den3$y,den4$y,den5$y))

pdf(file='samples_by100.pdf',width=5,height=5)

plot(den1$x,den1$y,
	# xlim=c(xmin,xmax),ylim=c(ymin,ymax),
	type='l',col='blue',xlab='1/L',ylab='probability density')
lines(den2$x,den2$y,col='red')
# lines(den3$x,den3$y,col='black')
# lines(den4$x,den4$y,col='green')
# lines(den5$x,den5$y,col='cyan')

empmode = den2$x[which.max(den2$y)]
abline(v=empmode)

abline(v=2*pi,lty=2)

# legend(x="topright",legend=c("h=by1","h=by2","h=by3"),col=c('blue','red','black','green', 'cyan'),lwd=1)

dev.off()
