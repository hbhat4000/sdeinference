diffdiffustheta1 <- function(c0, y, z)
{
	n = length(y)
	m = length(z)
	dtheta1g = rep(0,n)
	out = replicate(m,dtheta1g)
	return(out)
}