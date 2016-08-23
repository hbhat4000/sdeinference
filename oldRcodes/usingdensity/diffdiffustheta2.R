diffdiffustheta2 <- function(c0, y,z)
{
	n = length(y)
	m = length(z)
	dtheta2g = rep(1,n)
	out = replicate(m,dtheta2g)
	return(out)
}