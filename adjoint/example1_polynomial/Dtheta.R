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
    part1 = (X - Y - f*h) 

    # these derivatives are the same regardless of what f and g you have
    dGdf = G*part1/g^2
    dGdg = G*(-1/g + (part1)^2/(g^3*h))

    # make everything a list so that we can handle high-dimensional problems
    derivatives = list(NULL)
    dfdtheta = list(NULL)
    
    # chain rule!
    nc0 = length(c0)
    for (i in c(1:length(c0)))
   {

	  if (i < nc0)
        {
            dfdtheta =  Y^(i-1)
            dgdtheta = 0
        }
        else
        {
            dfdtheta = 0
            dgdtheta = 1
        }
	 
        derivatives[[i]] = dGdf * dfdtheta + dGdg * dgdtheta
   }
    
    return(derivatives)
}

