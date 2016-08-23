allsol = matrix(0,nrow=9,ncol=6)
allsol[1,] = c( 0.361, 0.968, 0.836, 39, 0.05, 50 )
allsol[2,] = c( 0.362, 0.968, 0.839, 46, 0.02, 50 )
allsol[3,] = c( 0.362, 0.968, 0.840, 42, 0.01, 50 )
allsol[4,] = c( 0.362, 0.968, 0.841, 28, 0.005, 50 )
allsol[5,] = c( 0.463, 0.885, 0.966, 45, 0.05, 300 )
allsol[6,] = c( 0.466, 0.886, 0.973, 22, 0.02, 300 )
allsol[7,] = c( 0.467, 0.886, 0.975, 22, 0.01, 300 )
allsol[8,] = c( 0.468, 0.886, 0.976, 26, 0.005, 300 )
allsol[9,] = c( 0.468, 0.886, 0.976, 20, 0.002, 300 )
truesol = c(0.5, 0.9, 1)

rmserror = numeric(length = nrow(allsol))
for (i in c(1:nrow(allsol)))
{
    rmserror[i] = sqrt(sum( (truesol - allsol[i,1:3])^2 ))
}

allsol = cbind(allsol, rmserror)


