## Application to real data
## CKLS model vs CIR model
## CKLS (mod1): dX(t) = (theta1+theta2* X(t))* dt + theta3 * X(t)^theta4 * dW(t)
## CIR (mod2): dX(t) = (theta1+theta2* X(t))* dt + theta3 * sqrt(X(t)) * dW(t)

set.seed(1234)
data(Irates)

rates = Irates[,"r1"]
rates = window(rates, start = 1964.471, end = 1989.333)

fx1 = expression(theta[1] + theta[2]*x)
gx1 = expression(theta[3]*x^theta[4])
gx2 = expression(theta[3]*sqrt(x))

init = c(1, 1, 1, 1)
fitmod1 = fitsde(rates, drift = fx1, diffusion = gx1, pmle = "euler", start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3], theta4 = init[4]), optim.method = "L-BFGS-B")
fitmod2 = fitsde(rates, drift = fx1, diffusion = gx2, pmle = "euler", start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), optim.method = "L-BFGS-B")

# summary(fitmod1)
# summary(fitmod2)
coef(fitmod1)
coef(fitmod2)
# confint(fitmod1, parm = c('theta2', 'theta3'))
# confint(fitmod2, parm = c('theta2', 'theta3'))
# AIC(fitmod1)
# AIC(fitmod2)

## Display
## CKLS Modele
op = par(mfrow = c(1, 2))
theta = coef(fitmod1)
N = length(rates)

# M: number of sample paths
res = snssde1d(drift = fx1, diffusion = gx1, M = 50, t0 = time(rates)[1], T = time(rates)[N], Dt = deltat(rates), x0 = rates[1], N)
plot(res, plot.type = "single", ylim = c(0,40))
lines(rates, col = 2, lwd = 2)
legend("topleft", c("real data", "CKLS model"), inset = .01, col = c(2, 1), lwd = 2, cex = 0.8)

## CIR Model
theta = coef(fitmod2)
res = snssde1d(drift = fx1, diffusion = gx2, M = 50, t0 = time(rates)[1], T = time(rates)[N], Dt = deltat(rates), x0 = rates[1], N)

plot(res, plot.type = "single", ylim=c(0, 40))
lines(rates, col = 2, lwd = 2)
legend("topleft", c("real data", "CIR model"), inset = .01, col = c(2, 1), lwd = 2, cex = 0.8)
par(op)