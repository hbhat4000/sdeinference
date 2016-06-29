Dtheta <- function(xvec, yvec, h, theta)
{
  # each column same and replicated, rows have xvec
  X = replicate(length(yvec), xvec)
  # each row same, columns have yvec
  Y = replicate(length(xvec), yvec)
  if (length(yvec)>1) Y = t(Y)
  else Y = as.matrix(Y)

  f = drift(theta, Y)
  g = abs(diffusion(theta, Y))
  part1 = (X - Y - f*h) 
  part2 = 2*(g^2)*h
  thisdiff = abs(diffusion(theta, Y))
  G = exp(-(part1)^2/(part2))/(sqrt(pi*part2))
  
  # these derivatives are the same regardless of what f and g you have
  dGdf = G*part1/g^2
  dGdg = G*(-1/g + (part1)^2/(g^3*h))
  
  # derivatives with respect to theta
  # make everything a list so that we can handle high-dimensional problems
  derivatives = list(NULL)
  dfdtheta = list(NULL)
  dgdtheta = list(NULL)
  dfdtheta[[1]] = f/theta[1]
  dgdtheta[[1]] = 0
  dfdtheta[[2]] = theta[1]
  dgdtheta[[2]] = 0
  dfdtheta[[3]] = 0
  dgdtheta[[3]] = 1
  
  # chain rule!
  for (i in c(1:3))
    derivatives[[i]] = dGdf * dfdtheta[[i]] + dGdg * dgdtheta[[i]]
    
    return(derivatives)
}


