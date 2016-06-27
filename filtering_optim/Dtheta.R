Dtheta <- function(xvec, yvec, h, driftfun, difffun1, c0)
{
  nn = length(xvec)
  mm = length(yvec)
  
  # each column same and replicated, rows have xvec
  X = replicate(length(yvec), xvec)
  # each row same, columns have yvec
  Y = t(replicate(length(xvec), yvec))
  
  f = driftfun(c0 ,Y)
  g = abs(difffun1(c0, Y))
  part1 = (X - Y - f*h) 
  part2 = 2*(g^2)*h
  thisdiff = abs(difffun1(c0, Y))
  G = exp(-(part1)^2/(part2))/(sqrt(pi*part2))
  
  # G = integrandmat(xvec, yvec, h, driftfun, difffun1, c0)
  # X = replicate(mm, xvec)
  # Y = t(replicate(nn, yvec))
  
  # these derivatives are the same regardless of what f and g you have
  dGdf = G*part1/g^2
  dGdg = G*(-1/g + (part1)^2/(g^3*h))
  
  # derivatives with respect to theta
  # make everything a list so that we can handle high-dimensional problems
  derivatives = list(NULL)
  dfdtheta = list(NULL)
  dgdtheta = list(NULL)
  dfdtheta[[1]] = f/c0[1]
  dgdtheta[[1]] = 0
  dfdtheta[[2]] = c0[1]
  dgdtheta[[2]] = 0
  dfdtheta[[3]] = 0
  dgdtheta[[3]] = 1
  
  # chain rule!
  for (i in c(1:3))
    derivatives[[i]] = dGdf * dfdtheta[[i]] + dGdg * dgdtheta[[i]]
  
  # derivatives with respect to x
  derivatives[[4]] = -2*(part1)*G/(part2)
  
  # derivatives with respect to sigeps2
  derivatives[[5]] = 
    
    return(derivatives)
}


