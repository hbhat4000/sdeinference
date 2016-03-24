Dtheta <- function(xvec, yvec, h, driftfun, difffun1, c0)
{
	G = integrandmat(xvec, yvec, h, driftfun, difffun1, c0)
	nn = length(xvec)
    mm = length(yvec)
  
    # each column same and replicated, rows have xvec
	X = replicate(mm, xvec)
	
	# each row same, columns have yvec
 	Y = t(replicate(nn, yvec))
 	
 	f = driftfun(c0 ,Y)
	g = difffun1(c0, Y)
	
 	part1 = (X - Y - f*h)    # (x - mu) on grid
	part2 = g^2
	
	derivative1 = G * (part1) * (f / c0[1]) / (part2)
	derivative2 = G * (part1) * (c0[1] * Y) / (part2)
	derivative3 = G * ((part1)^2 - h*part2) / (h*part2*g) 
	
 	return(list(derivative1, derivative2, derivative3))
}
