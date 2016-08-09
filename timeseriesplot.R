xtraj = c(1:25)
sinval = matrix(0, nrow = 100, ncol = 25)
for(i in c(1:100))
{
  sinval[i,] = sin(xtraj) + rnorm(n = length(xtraj), 0, 0.1)
}

plot(sinval[1,], type = "l", col = "azure4")
for(i in c(2:100))
{
  lines(sinval[i,], col = "azure4")
}