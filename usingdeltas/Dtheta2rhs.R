Dtheta2rhs <- function(xvec,h,driftfun,difffun1,c0,data)
{
	G = integrandmat(xvec,data,h,driftfun,difffun1,c0)
 	# diffmat = diffgandf(c0,xvec)
	nn = length(xvec)
	mm = length(data)
	X=replicate(mm,xvec)
 	Y=t(replicate(nn,data))
	internpart1 = (X-Y-driftfun(c0,Y)*h)
	# internpart2 = replicate(mm,diffmat[,4])
	internpart2 = diffdiffustheta2(c0,xvec, data)	
	part1 = -internpart2/difffun1(c0,Y)
	# part2 = difffun1(c0,Y)*internpart1*replicate(mm,diffmat[,2])
	part2 = difffun1(c0,Y)*internpart1*diffdrifttheta2(c0,xvec, data)
	part3 = internpart1^2*internpart2
	part4 = difffun1(c0,Y)^3*h
	out = G*(part1+(part2+part3)/part4)
 	return(out)
}
