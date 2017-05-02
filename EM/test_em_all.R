rm(list = ls(all = TRUE))

library("matrixcalc")

# xtraj = matrix(nrow = 1, ncol = 2)

# load necessary functions
source('dtq_main.R')
source('Dtheta.R')
# source('kolmogorov_compare.R')

k = 0.01
M = 800
deltat = 1
numsteps = 50
h = deltat/numsteps
theta = c(1, 0, 2)
init = 1
final = 3

allout = dtq_all(theta, h, k, M, numsteps, init, final)
print(c(sum(allout$first)*k,sum(allout$last)*k))

# exactcomplete = transition_forward(theta,x=final,y=init,t=deltat)
# print(c(log(allout$complete),log(exactcomplete)))

# part1 = transition_forward(theta,x=grid,y=init,t=h)
# part2 = transition_forward(theta,x=final,y=grid,t=(deltat-h))
# exactfirst = hadamard.prod(part1,part2)/exactcomplete

# print(k*sum(abs(allout$first - exactfirst)))

# part1 = transition_forward(theta,x=grid,y=init,t=(deltat-h))
# part2 = transition_forward(theta,x=final,y=grid,t=h)
# exactlast = hadamard.prod(part1,part2)/exactcomplete

# print(k*sum(abs(allout$last - exactlast)))

# gridmat = replicate(length(grid), grid)
# part2 = transition_forward(theta,x=gridmat,y=t(gridmat),t=h)
# 
# for (j in c(1:(numsteps-2)))
# {
#   print(sum(allout$pdf2d[[j]])*k^2)
#   part1 = matrix(transition_forward(theta,x=grid,y=init,t=(j*h)),ncol=1)
#   part3 = matrix(transition_forward(theta,x=final,y=grid,t=(deltat - (j+1)*h)),nrow=1)
#   exactinternal = hadamard.prod(t(kronecker(part3,part1)),part2) / exactcomplete
#   print(sum(exactinternal)*k^2)
#   print(k^2*sum(abs(allout$pdf2d[[j]] - exactinternal)))
# }
