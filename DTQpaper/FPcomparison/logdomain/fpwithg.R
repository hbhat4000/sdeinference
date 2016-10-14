fperr <- function(dt,T,s,driftfun,difffun,exactfun)
{
  ptm = proc.time()
  dx=dt^s
  L=(s+2)*(-log(dt))
  xgrid=seq(from=-L,to=(L+dx),by=dx)
  m=length(xgrid)-3
  x=xgrid[2:(length(xgrid)-1)]
  m=length(x)-1
  nu=dt/dx
  
  Nt=ceiling(T/dt)
  v=numeric(length=length(x))
  
  f = driftfun(xgrid)
  g = difffun(xgrid)
  g2=g^2
  kappa=1 # used to be mean(g2) but that leads to divergent code
  
  r1=-0.5*nu*f[3:(m+2)]
  l1=0.5*nu*f[2:(m+1)]
  Bsparse=bandSparse(m+1,k=c(-1,0,1),diagonals=list(l1,rep(1,m+1),r1))
  
  nu2=dt/(dx^2)
  c2=rep(1,m+1)+nu2*g2[2:(m+2)]
  r2=-0.5*nu2*g2[3:(m+2)]
  l2=-0.5*nu2*g2[2:(m+1)]
  Asparse=bandSparse(m+1,k=c(-1,0,1),diagonals=list(l2,c2,r2))
  
  C = solve(Asparse,Bsparse)
  
  #print(max(abs(eigen(C)$values)))
  
  for (n in c(1:Nt))
  {
    t=(n-1)*dt
    u=(1/sqrt(2*pi*kappa*n*dt))*exp(-(xgrid^2)/(2*kappa*n*dt))
    gvec1=0.5*nu*(f[1:(m+1)]*u[1:(m+1)] - f[3:(m+3)]*u[3:(m+3)])
    ga=(g2-kappa)*u
    gvec2=0.5*nu2*(ga[1:(m+1)] - 2*ga[2:(m+2)] + ga[3:(m+3)])
    v = C%*%v + solve(Asparse,(gvec1+gvec2))
    u=u[2:(m+2)]
    p=u+v
  }
  timetaken = proc.time() - ptm
  
  pexact = exactfun(x,T)
  l1err=sum(abs(p-pexact)*dx)
  return(list(timing=timetaken,error=l1err))
}



