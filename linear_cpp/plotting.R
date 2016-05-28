rm(list=ls(all=TRUE))
library('latex2exp')
hmcdatanum = 10
fname = paste('HMCsamples_', hmcdatanum, '.RData', sep = '')
load(fname)
attach(myout)

# acceptance = 1 - mean(duplicated(thetamat[-(1:burnin),]))
newartrack = artrack[-(1:burnin)]
acceptance = sum(artrack)/length(artrack)

########## First plotting technique ##########
# par(mfrow = c(2,3))
par(mfrow = c(2,2))

# hist(thetamat[-(1:burnin),1], nclass = 30, main = "Posterior of theta1", xlab = TeX("$\\theta_1$ values"), xlim = c(actualtheta[1] - 0.05, actualtheta[1] + 0.05))
# abline(v = mean(thetamat[-(1:burnin),1]), col = "green")
# abline(v = actualtheta[1], col = "red")
# lines(density(thetamat[-(1:burnin),2]), col = "blue")

hist(thetamat[-(1:burnin),2], nclass = 30, main = TeX("Posterior of $\\theta_2$"), xlab = TeX("$\\theta_2$ values"), xlim = c(actualtheta[2] - 0.05, actualtheta[2] + 0.05))
abline(v = mean(thetamat[-(1:burnin),2]), col = "green")
abline(v = actualtheta[2], col = "red")
lines(density(thetamat[-(1:burnin),2]), col = "blue")

# ## For theta_3^2 as the parameter ###
# thetamat3new = (thetamat[-(1:burnin), 3])^2
# newtheta3 = (actualtheta[3])^2
# hist(thetamat3new, nclass = 30, main = TeX("Posterior of $\\theta_3^2$"), xlab = TeX("$\\theta_3^2$ values"), xlim = c(newtheta3 - 0.05, newtheta3 + 0.05))
# abline(v = mean(thetamat3new), col = "green")
# abline(v = newtheta3, col = "red")
# lines(density(thetamat3new), col = "blue")

hist(thetamat[-(1:burnin),3], nclass = 30, main = TeX("Posterior of $\\theta_3$"), xlab = TeX("$\\theta_3$ values"), xlim = c(actualtheta[3] - 0.05, actualtheta[3] + 0.05))
abline(v = mean(thetamat[-(1:burnin),3]), col = "green")
abline(v = actualtheta[3], col = "red")
lines(density(thetamat[-(1:burnin),3]), col = "blue")

# plot(thetamat[-(1:burnin),1], type = "l", xlab = "Iterations", ylab = TeX("$\\theta_1$ values"), main = TeX("Trace plot of $\\theta_1$"), ylim = c(actualtheta[1] - 0.05, actualtheta[1] + 0.05))
# abline(h = actualtheta[1], col = "red")

plot(thetamat[-(1:burnin),2], type = "l", xlab = "Iterations", ylab = TeX("$\\theta_2$ values"), main = TeX("Trace plot of $\\theta_2$"), ylim = c(actualtheta[2] - 0.05, actualtheta[2] + 0.05))
abline(h = actualtheta[2], col = "red")

# ## For theta_3^2 as the parameter ###
# plot(thetamat3new, type = "l", xlab = "Iterations", ylab = TeX("$\\theta_3^2$ values"), main = TeX("Trace plot of $\\theta_3^2$"), ylim = c(newtheta3 - 0.05, newtheta3 + 0.05))
# abline(h = newtheta3, col = "red")

plot(thetamat[-(1:burnin),3], type = "l", xlab = "Iterations", ylab = TeX("$\\theta_3$ values"), main = TeX("Trace plot of $\\theta_3$"), ylim = c(actualtheta[3] - 0.05, actualtheta[3] + 0.05))
abline(h = actualtheta[3], col = "red")

# for comparison:
# summary(lm(y~x))

########## Second plotting technique ##########
# Using MCMC package for plotting
# library('coda')
# chain = mcmc(thetamat)
# summary(chain)
# plot(chain)
