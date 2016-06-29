Dtheta <- function(xvec, yvec, h, theta)
{
  # each column same and replicated, rows have xvec
  X = replicate(length(yvec), xvec)
  # each row same, columns have yvec
  Y = t(replicate(length(xvec), yvec))
  
  f = drift(theta, Y)
  g = abs(diffusion(theta, Y))
  part1 = (X - Y - f*h) 
  part2 = 2*(g^2)*h
  G = exp(-(part1)^2/(part2))/(sqrt(pi*part2))
    
  return(G)
}


