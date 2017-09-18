# Note that this code is a modified version of the code that appears on
# Dr. Yoon Dong Lee's web site: http://blog.naver.com/widylee/120128541187
rm(list=ls(all=TRUE))

# Harish S. Bhat's bistable SDE
# we have included the Lamperti transformation to unit diffusion coefficient
hsb.model <- list(mux=quote(p1*x*(p2-x^2)), sigmax=quote(p3*sqrt(1+p4^2*x^2)), trf=quote(asinh(p4*x)/(p3*p4)))

# The working environment w.e is used to pass the arguments in evaluating symbloically
# represented functions.
# my.Dmuy will contain the symbolic functions of mu.yox and its derivatives.
# my.Model will contain the model representation.
# my.Pams will contains the string of the paramer names used in the model.
# my.Mname will contain the model names, such as u0, u1, cir, etc.
# my.KK will contains the maximum allowed K values.

w.e<-new.env()
my.Dmuy<- NULL
my.Model<-NULL
my.Pams<-NULL
my.Mname<-NULL
my.KK<-4

# The function set set.model sets the variables my.Dmuy and my.Pams, etc.
set.model<-function(model) {
  
  mu.yox<-function(){ 
    zx<-quote(a/b-c/2)
    zx[[2]][[2]]<- model$mux
    zx[[2]][[3]]<- model$sigmax
    zx[[3]][[2]]<- D(model$sigmax,"x")
    zx
  }
  
  # mu.yox gives the symbolic function of mu of Y_t as a function of x
  # Dy gives the symbolic derivatives w.r.t. y  of a function of x,
  
  Dy<-function(muyx){
    zx<-quote(a*b)
    zx[[2]]<- D(muyx,"x")     
    zx[[3]]<- model$sigmax
    zx
  }
  
  my.Mname<<-substitute(model)
  my.Model<<- model
  
  # The following steps extract names of parameter used in the model description
  # All parameters should have the form "p+numbers" such as "p0", "p1", "p01" and "p21".
  
  tmp<-all.vars(model$mux)
  tmp<-c(tmp, all.vars(model$sigmax) )
  my.Pams<<-unique(sort(tmp[grep("^p[0-9]+$",tmp)]))
  
  # my.Dmuy will contain symbolic derivatives of mu.yox upto (2*KK-1) order
  
  my.Dmuy<<-list(mu.yox())
  for(i in 1:(2*(my.KK-1))) my.Dmuy[[length(my.Dmuy)+1]]<<-Dy(my.Dmuy[[length(my.Dmuy)]])
  
}



# You need to specify your model, here, by directly writing the function mux and sigmax.
# Otherwise, you can specify the model by using prespecified models as follows.


likelihood<- function( p, xx, delta, J=4, K=3 ){
  
  # Hermite polynomials H_j(z) and  H_j(0)
  
  hermite <- function(j,x){
    if(j<0) res<-rep(0, length(x)) 
    if(j==0) res<-rep(1,length(x)) 
    if(j==1) res<-x
    if (j>1) res <- x*hermite(j-1,x)-(j-1)*hermite(j-2,x)
    res
  }
  
  hermite0 <- function(j) {
    res<-0
    m<- trunc(j/2)
    if(m>=0){if( m*2!=j ) res<-0  else res<-prod(1-2*(0:m))}
    res
  }  
  
  ut7<- matrix( c(
    3,2,0,0,0,0,0, 1,3,0,0,0,0,0, 4,0,1,0,0,0,0, 2,1,1,0,0,0,0, 0,2,1,0,0,0,0,
    1,0,2,0,0,0,0, 3,0,0,1,0,0,0, 1,1,0,1,0,0,0, 0,0,1,1,0,0,0, 2,0,0,0,1,0,0,
    0,1,0,0,1,0,0, 1,0,0,0,0,1,0, 0,0,0,0,0,0,1), 7, )
  
  ut6<- matrix( c(
    4,1,0,0,0,0, 2,2,0,0,0,0, 0,3,0,0,0,0, 3,0,1,0,0,0, 1,1,1,0,0,0, 0,0,2,0,0,0,
    2,0,0,1,0,0, 0,1,0,1,0,0, 1,0,0,0,1,0, 0,0,0,0,0,1), 6, )
  
  ut5<- matrix( c(
    5,0,0,0,0, 3,1,0,0,0, 1,2,0,0,0, 2,0,1,0,0, 0,1,1,0,0, 1,0,0,1,0, 0,0,0,0,1),5,)
  
  ut4<- matrix( c(4,0,0,0, 2,1,0,0, 0,2,0,0, 1,0,1,0, 0,0,0,1), 4,)
  
  ut3<- matrix( c( 3,0,0, 1,1,0, 0,0,1), 3,)
  
  ut2<- matrix( c(2,0,0,1),2,)
  
  ct7<-c(0,1,0,4,7/2,7/2,1,11/2,11/4,3/2,7/4,3/4,1/8)
  ct6<-c(0,7,4,4,22,21/4,7,8,4,3/4)
  
  ct5<-function(k) c( (k-3)*(k-4)/60, (k-3)*(2*k-5)/12, (k-2)*(3*k-8)/12, (4*k^2 -17*k+17)/24,
                      (2*k-3)*(2*k-5)/24, (k-1)*(2*k-5)/24, (4*k^2 -13*k +8)/240 )*k*(k-1)*(k-2)/(2^(k-4))
  
  ct4<-function(k) c( (k-3)/12, (k-2)/2, (3*k-5)/12, (4*k-5)/12, (k-1)/12 )*k*(k-1)*(k-2)/(2^(k-3))
  ct3<-function(k) c( (k-2)/3 , (2*k-3)/2, (4*k-5)/12 )*k*(k-1)/(2^(k-2))
  ct2<-function(k) c(k*(k-1)/(2^(k-1)), k*(k-1)/(2^(k-1)))
  
  cts<-function(k,i) {
    if(2*k-i==0) res<-1/(2^k)
    if(2*k-i==1) res<-k/(2^(k-1))
    if(2*k-i==2) res<-ct2(k)
    if(2*k-i==3) res<-ct3(k)
    if(2*k-i==4) res<-ct4(k)
    if(2*k-i==5) res<-ct5(k)
    if(2*k-i==6) res<-ct6
    if(2*k-i==7) res<-ct7
    res
  }
  
  uts<-function(k,i) { 
    if(2*k-i==0) res<-matrix(0,1,1) 
    if(2*k-i==1) res<-matrix(1,1,1)
    if(2*k-i==2) res<-ut2 
    if(2*k-i==3) res<-ut3
    if(2*k-i==4) res<-ut4
    if(2*k-i==5) res<-ut5
    if(2*k-i==6) res<-ut6
    if(2*k-i==7) res<-ut7
    res
  }
  
  ## uu(k,i) in this program is the uu(k-1,i)
  ## in the paper Lee, et.al. (2010)
  
  uu<-function(k,i){
    if(2*k-i ==0) res<-cts(k,i)  
    if(2*k-i ==1) res<-cts(k,i)*eval(my.Dmuy[[1]],w.e)
    if(2*k-i >1){
      tu<-uts(k,i)
      tc<-cts(k,i)
      res<-0
      for(l in 1:length(tc)){ 
        dord<-(1:(2*k-i))[tu[,l]>0]
        pord<-tu[tu[,l]>0,l]
        pres<-tc[l]    
        for(r in 1:length(dord)){ 
          if(all(pres==0)) break ;
          pres<-pres*(eval(my.Dmuy[[dord[r]]], w.e)**pord[r])  }               
        res<-res+pres
      }
    }
    res
  }
  
  
  etaz<-function(j,K,delta){
    res<-hermite0(j)/factorial(j)
    if(K>=1) for(k in 1:K) for(i in 1:(2*k)) if(i<=j) res<-res+
        uu(k,i)*hermite0(j-i)*(delta^(k-i/2))/(factorial(k)*factorial(j-i))
    res
  }
  
  
  hermite.series<-function(z,J,K,delta){
    res<-0
    for(j in 0:J) res<-res+hermite(j,z)*etaz(j,K,delta)
    res
  }
  
  
  # recip.sigmax will be used for numerical integration when "trf" is not given  
  
  recip.sigmax<-function(v){
    assign("x", v, w.e)
    1 / eval(my.Model$sigmax, w.e)
  }
  
  
  K<- min(K,my.KK) # for preventing K from going over my.KK
  
  if(length(p)!=length(my.Pams)) error<-1
  for(i in 1:min(length(p), length(my.Pams))) assign(my.Pams[i], p[i], w.e )
  
  
  x0<-xx[-length(xx)]
  x1<-xx[-1]  
  
  assign( "x" , x1,  w.e)
  y1<-eval(my.Model$trf,  w.e)
  sigma1<-eval(my.Model$sigmax,  w.e)
  
  # When "trf" is specified, evaluation from symbolic representaion,
  # otherwise numerical integration will be used.
  
  if(!is.null(my.Model$trf)){assign( "x" , x0,  w.e); y0<-eval(my.Model$trf, w.e); z<-(y1-y0)/sqrt(delta) }
  else { z<-rep(0,length(x1))
  for(i in 1:length(x1)) z[i]<- integrate(recip.sigmax ,x0[i], x1[i])$value
  z<- z/sqrt(delta)
  assign( "x" , x0,  w.e)
  }
  
  dnorm(z)*hermite.series(z,J,K,delta)/(sigma1*sqrt(delta)) 
  
}

batchlikelihood<- function( pmat, xx, delta, J=4, K=3 ){
  
  # Hermite polynomials H_j(z) and  H_j(0)
  
  hermite <- function(j,x){
    if(j<0) res<-rep(0, length(x)) 
    if(j==0) res<-rep(1,length(x)) 
    if(j==1) res<-x
    if (j>1) res <- x*hermite(j-1,x)-(j-1)*hermite(j-2,x)
    res
  }
  
  hermite0 <- function(j) {
    res<-0
    m<- trunc(j/2)
    if(m>=0){if( m*2!=j ) res<-0  else res<-prod(1-2*(0:m))}
    res
  }  
  
  ut7<- matrix( c(
    3,2,0,0,0,0,0, 1,3,0,0,0,0,0, 4,0,1,0,0,0,0, 2,1,1,0,0,0,0, 0,2,1,0,0,0,0,
    1,0,2,0,0,0,0, 3,0,0,1,0,0,0, 1,1,0,1,0,0,0, 0,0,1,1,0,0,0, 2,0,0,0,1,0,0,
    0,1,0,0,1,0,0, 1,0,0,0,0,1,0, 0,0,0,0,0,0,1), 7, )
  
  ut6<- matrix( c(
    4,1,0,0,0,0, 2,2,0,0,0,0, 0,3,0,0,0,0, 3,0,1,0,0,0, 1,1,1,0,0,0, 0,0,2,0,0,0,
    2,0,0,1,0,0, 0,1,0,1,0,0, 1,0,0,0,1,0, 0,0,0,0,0,1), 6, )
  
  ut5<- matrix( c(
    5,0,0,0,0, 3,1,0,0,0, 1,2,0,0,0, 2,0,1,0,0, 0,1,1,0,0, 1,0,0,1,0, 0,0,0,0,1),5,)
  
  ut4<- matrix( c(4,0,0,0, 2,1,0,0, 0,2,0,0, 1,0,1,0, 0,0,0,1), 4,)
  
  ut3<- matrix( c( 3,0,0, 1,1,0, 0,0,1), 3,)
  
  ut2<- matrix( c(2,0,0,1),2,)
  
  ct7<-c(0,1,0,4,7/2,7/2,1,11/2,11/4,3/2,7/4,3/4,1/8)
  ct6<-c(0,7,4,4,22,21/4,7,8,4,3/4)
  
  ct5<-function(k) c( (k-3)*(k-4)/60, (k-3)*(2*k-5)/12, (k-2)*(3*k-8)/12, (4*k^2 -17*k+17)/24,
                      (2*k-3)*(2*k-5)/24, (k-1)*(2*k-5)/24, (4*k^2 -13*k +8)/240 )*k*(k-1)*(k-2)/(2^(k-4))
  
  ct4<-function(k) c( (k-3)/12, (k-2)/2, (3*k-5)/12, (4*k-5)/12, (k-1)/12 )*k*(k-1)*(k-2)/(2^(k-3))
  ct3<-function(k) c( (k-2)/3 , (2*k-3)/2, (4*k-5)/12 )*k*(k-1)/(2^(k-2))
  ct2<-function(k) c(k*(k-1)/(2^(k-1)), k*(k-1)/(2^(k-1)))
  
  cts<-function(k,i) {
    if(2*k-i==0) res<-1/(2^k)
    if(2*k-i==1) res<-k/(2^(k-1))
    if(2*k-i==2) res<-ct2(k)
    if(2*k-i==3) res<-ct3(k)
    if(2*k-i==4) res<-ct4(k)
    if(2*k-i==5) res<-ct5(k)
    if(2*k-i==6) res<-ct6
    if(2*k-i==7) res<-ct7
    res
  }
  
  uts<-function(k,i) { 
    if(2*k-i==0) res<-matrix(0,1,1) 
    if(2*k-i==1) res<-matrix(1,1,1)
    if(2*k-i==2) res<-ut2 
    if(2*k-i==3) res<-ut3
    if(2*k-i==4) res<-ut4
    if(2*k-i==5) res<-ut5
    if(2*k-i==6) res<-ut6
    if(2*k-i==7) res<-ut7
    res
  }
  
  ## uu(k,i) in this program is the uu(k-1,i)
  ## in the paper Lee, et.al. (2010)
  
  uu<-function(k,i){
    if(2*k-i ==0) res<-cts(k,i)  
    if(2*k-i ==1) res<-cts(k,i)*eval(my.Dmuy[[1]],w.e)
    if(2*k-i >1){
      tu<-uts(k,i)
      tc<-cts(k,i)
      res<-0
      for(l in 1:length(tc)){ 
        dord<-(1:(2*k-i))[tu[,l]>0]
        pord<-tu[tu[,l]>0,l]
        pres<-tc[l]    
        for(r in 1:length(dord)){ 
          if(all(pres==0)) break ;
          pres<-pres*(eval(my.Dmuy[[dord[r]]], w.e)**pord[r])  }               
        res<-res+pres
      }
    }
    res
  }
  
  uut = array(0, dim=c(K,2*K,(length(xx)-1)))
  
  etaz<-function(j,K,delta){
    res<-hermite0(j)/factorial(j)
    if(K>=1) for(k in 1:K) for(i in 1:(2*k)) if(i<=j) res<-res+
        uut[k,i,]*hermite0(j-i)*(delta^(k-i/2))/(factorial(k)*factorial(j-i))
    res
  }
  
  
  hermite.series<-function(z,J,K,delta){
    res<-0
    for(j in 0:J) res<-res+hermite(j,z)*etaz(j,K,delta)
    res
  }
  
  
  # recip.sigmax will be used for numerical integration when "trf" is not given  
  
  recip.sigmax<-function(v){
    assign("x", v, w.e)
    1 / eval(my.Model$sigmax, w.e)
  }
  
  
  K<- min(K,my.KK) # for preventing K from going over my.KK
  
  outmat = matrix(0, nrow = nrow(pmat),ncol = (length(xx)-1))
  x0<-xx[-length(xx)]
  x1<-xx[-1]
  for (iii in c(1:nrow(pmat)))
  {
    p = pmat[iii,]
    if(length(p)!=length(my.Pams)) error<-1
    for(i in 1:min(length(p), length(my.Pams))) assign(my.Pams[i], p[i], w.e )
    
    assign( "x" , x1,  w.e)
    y1<-eval(my.Model$trf,  w.e)
    sigma1<-eval(my.Model$sigmax,  w.e)

    # When "trf" is specified, evaluation from symbolic representaion,
    # otherwise numerical integration will be used.
    
    if(!is.null(my.Model$trf)){
      assign( "x" , x0,  w.e)
      y0<-eval(my.Model$trf, w.e)
      z<-(y1-y0)/sqrt(delta)
    }
    else {
      z<-rep(0,length(x1))
      for(i in 1:length(x1)) z[i]<- integrate(recip.sigmax ,x0[i], x1[i])$value
      z<- z/sqrt(delta)
      assign( "x" , x0,  w.e)
    }
    for (k in 1:K)
      for (i in 1:(2*k))
        uut[k,i,] = uu(k,i)
    outmat[iii,] = dnorm(z)*hermite.series(z,J,K,delta)/(sigma1*sqrt(delta))
  }
  outmat
}

#### ###### ###### ####### ############ ########## ########## ###########

####### ############ ########## ########## ########### ############# 


# cir model has exact likelihood function. likelihood.cir gives the the values.

likelihood.cir<-function( p, xx, delta){
  x0<-xx[-length(xx)]; x1<-xx[-1]
  alpha<-p[1]; beta<- p[2]; sigma<- p[3]
  alsig2<-  alpha/(sigma^2)
  aldelta<- alpha*delta
  cz<- 2* alsig2/(1-exp(-aldelta))
  x0tilde<- x0*exp(-aldelta)
  kk<- 2*alsig2*beta-1
  tmp<- cz*((x1/x0tilde)**(kk/2))*exp(-cz*((sqrt(x1)-sqrt(x0tilde))**2))
  tmp*besselI(2*cz*sqrt(x1*x0tilde),kk,exp=TRUE)
}


# For generating diffusion process, burning step is necessary.

burning<-function(x, n , delta ){
  assign("x",x, w.e)
  for(i in 1:n) assign("x",
                       get("x",w.e) + eval(my.Model$mux, w.e)*delta+ eval(my.Model$sigmax,w.e)*sqrt(delta)*rnorm(1),
                       w.e)
  x
}

# For generating diffusion process, my.Model and my.Pams are saved in set.param

set.param<-function(p,model){
  my.Model<<- model
  tmp<-all.vars(model$mux)
  tmp<-c(tmp, all.vars(model$sigmax) )
  my.Pams<<-unique(sort(tmp[grep("^p[0-9]+$",tmp)]))  
  if(length(p)!=length(my.Pams)) error<-1
  for(i in 1:min(length(p),length(my.Pams))) assign(my.Pams[i], p[i], w.e)      
}

# n1 is the size (length) of sample to be generated
# n2 is the number of subintervals dividing an interval of observations.

gen<-function(x.init, delta, n1=size.sample, n2=step.gen){
  x.gen<-NULL
  # burning(burning(x.init, 100, delta), 200 , delta/n2 )
  assign("x",x.init, w.e)
  for(i in 1:n1){  for(j in 1:n2)  assign("x",
                                          x<- get("x", w.e) + eval(my.Model$mux,w.e)*(delta/n2) + eval(my.Model$sigmax,w.e)*sqrt(delta/n2)*rnorm(1),
                                          w.e )
    x.gen<-c(x.gen,x) }
  x.gen
}



############# ######## ######### ############ ######### ############

# ACTUAL CODE THAT SWEEPS OUT LIKELIHOOD SURFACE

pp<-c(1,8,0.25,3.0)
delta<-0.1

set.model(hsb.model)
n0 = 50
n1 = 50
theta0grid = seq(from=-1,to=4,length.out=n0)
theta1grid = seq(from=5,to=12.5,length.out=n1)

ppmat = matrix(0,nrow=n0*n1,ncol=length(pp))
counter = 1
for (i0 in c(1:n0))
{
  for (i1 in c(1:n1))
  {
    pp[1] = theta0grid[i0]
    pp[2] = theta1grid[i1]
    ppmat[counter,] = pp
    counter = counter + 1
  }
}

Jvec = c(4, 8, 24)
Kvec = c(3, 4, 4)

xr1 = as.numeric(read.table('../traj1.txt',sep=',',header=FALSE))
xr2 = as.numeric(read.table('../traj2.txt',sep=',',header=FALSE))

for (i in c(1:3))
{
    lltest = rowSums(batchlikelihood(ppmat, xr1, delta, J=Jvec[i], K=Kvec[i]))
    lltest = lltest + rowSums(batchlikelihood(ppmat, xr2, delta, J=Jvec[i], K=Kvec[i]))
    fname = paste('llcf_',Jvec[i],'_',Kvec[i],'.RData',sep='')
    save.image(file=fname)
}

