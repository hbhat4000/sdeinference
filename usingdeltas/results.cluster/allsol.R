allsol = matrix(0,nrow=4,ncol=6)
allsol[1,] = c( 0.361, 0.968, 0.836, 39, 0.05, 50 )
allsol[2,] = c( 0.362, 0.968, 0.839, 46, 0.02, 50 )
allsol[3,] = c( 0.362, 0.968, 0.840, 42, 0.01, 50 )
allsol[4,] = c( 0.362, 0.968, 0.841, 28, 0.005, 50 )
truesol = c(0.5, 0.9, 1)

rmserror = numeric(length = nrow(allsol))
for (i in c(1:nrow(allsol)))
{
    rmserror[i] = sqrt(sum( (truesol - allsol[i,1:3])^2 ))
}


