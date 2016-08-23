diffgandf <- function(c0,y)
{
	
	n=length(y)
	mat = matrix(0,nrow=n, ncol=4)
	dtheta1f = rep(0.5,n)
	mat[,1] = dtheta1f
	dtheta2f = rep(0,n)
	mat[,2] = dtheta2f
	dtheta1g = rep(0,n)
	mat[,3] = dtheta1g
	dtheta2g = rep(1,n)
	mat[,4] = dtheta2g
	return(mat)

}

