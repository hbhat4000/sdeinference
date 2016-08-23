diffdrifttheta2 <- function(c0, y, z)
{
	n = length(y)
	m =length(z)
	dtheta2f = rep(0,n)
	out = replicate(m,dtheta2f)
	return(out)
}