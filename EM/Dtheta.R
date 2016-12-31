Dtheta <- function(x, y, h, theta)
{
  fval = f(theta, y)
  gval = abs(g(theta, y))
  part1 = (x - y - fval*h)
  G = exp(- part1^2/(2 * gval^2 * h))/(gval * sqrt(2*pi*h))
  
  # these derivatives are the same regardless of what f and g you have
  dGdf = G*part1/gval^2
  dGdg = G*(-1/gval + (part1)^2/(g^3*h))
  
  # derivatives with respect to theta
  # make everything a list so that we can handle high-dimensional problems
  derivatives = list(NULL)
  dfdtheta = list(NULL)
  dgdtheta = list(NULL)
  dfdtheta[[1]] = fval/theta[1]
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


