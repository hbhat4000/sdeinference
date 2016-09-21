rm(list = ls(all = TRUE))
# Tests with Ornstein Uhlenback model
# dX = theta[1]*(theta[2] - X) dt + theta[3] dWt

# Data Simulation with truetheta
truetheta = c(0.5, 0, 1)

f = expression(truetheta[1]*(truetheta[2] - x))
g = expression(exp(truetheta[3]))
sim = snssde1d(drift = f, diffusion = g, x0 = 1, Dt = 0.01)
mydata = sim$X

# Inference
fx = expression(theta[1]*(theta[2]- x))
gx = expression(exp(theta[3]))

init = c(1, 2, 0.5)

# fitsde(data, drift, diffusion, start = list(), pmle = c("euler","kessler", "ozaki", "shoji"), optim.method = "L-BFGS-B", lower = -Inf, upper = Inf, ...)
fitEuler = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "euler")
fitOzaki = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "ozaki")
fitShoji = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "shoji", lower = c(0, 0, 0), upper = c(10, 10, 5))
fitKessler = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "kessler")

# Results
print(coef(fitEuler))
print(coef(fitOzaki))
print(coef(fitShoji))
print(coef(fitKessler))

# summary(fitmod)
# confint(fitmod, parm = c("theta1", "theta2"), level = 0.95)

## 1-dim fitsde
data(Irates)
rates <- Irates[,"r1"]
rates <- window(rates, start=1964.471, end=1989.333)
fx <- expression(theta[1]+theta[2]*x)
gx <- expression(theta[3]*x^theta[4])
## theta = (theta1,theta2,theta3,theta4), p=4
fitmod <- fitsde(rates,drift=fx,diffusion=gx,pmle="euler",start = list(theta1=1,theta2=1,theta3=1,theta4=1),optim.method = "L-BFGS-B")
print(fitmod)
print(summary(fitmod))
print(coef(fitmod))
print(logLik(fitmod))
print(AIC(fitmod))
print(BIC(fitmod))
print(vcov(fitmod))
print(confint(fitmod,level=0.95))