# OU process: dX = -mu dt + sigma dWt
# Vasicek process: dX = mu(theta - X) dt + sigma dWt

rm(list = ls(all = TRUE))
# Tests with Ornstein Uhlenback model
# dX = theta[1]*(theta[2] - X) dt + theta[3] dWt

# Data Simulation with truetheta
truetheta = c(0.5, 0, 1)

f = expression(truetheta[1]*(truetheta[2] - x))
g = expression(exp(truetheta[3]))

# snssde1d(N = 1000, M = 1, x0 = 0, t0 = 0, T = 1, Dt, drift, diffusion, alpha = 0.5, mu = 0.5, type = c("ito", "str"), method = c("euler", "milstein", "predcorr", "smilstein", "taylor", "heun", "rk1", "rk2", "rk3"), ...)
sim = snssde1d(drift = f, diffusion = g, M = 100, x0 = 1, Dt = 0.01)
mydata100 = sim$X

# Inference
fx = expression(theta[1]*(theta[2]- x))
gx = expression(exp(theta[3]))

init = c(1, 2, 0.5)

# fitsde(data, drift, diffusion, start = list(), pmle = c("euler","kessler", "ozaki", "shoji"), optim.method = "L-BFGS-B", lower = -Inf, upper = Inf, ...)
fitEuler = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "euler")
fitOzaki = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "ozaki")
fitShoji = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "shoji")
fitKessler = fitsde(data = mydata, drift = fx, diffusion = gx, start = list(theta1 = init[1], theta2 = init[2], theta3 = init[3]), pmle = "kessler")

# Results
print(coef(fitEuler))
print(coef(fitOzaki))
print(coef(fitShoji))
print(coef(fitKessler))

# summary(fitmod)
# confint(fitmod, parm = c("theta1", "theta2"), level = 0.95)