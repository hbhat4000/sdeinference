rm(list=ls(all=TRUE))
library('pomp')
library('magrittr')
library('dplyr')

# for reproducibility
set.seed(1)

load('./vanderpol/fakedata_vanderpol_fullres.RData')
mydata = X[seq(from=1,to=dim(X),by=100),]
mymod.dat = data.frame(t=mydata[,3],Y1=as.numeric(mydata[,1]),Y2=as.numeric(mydata[,2]))

ptm = proc.time()

timeinc = mydata[2,3] - mydata[1,3]
myh = timeinc/1

step.fun <- Csnippet("
 double dW1 = rnorm(0,sqrt(dt));
 double dW2 = rnorm(0,sqrt(dt));
 X1 += theta1*theta1*X2*dt + exp(-X1*X1)*0.25*dW1;
 X2 += -theta2*theta2*X1*dt + theta3*theta3*(1-X1*X1)*X2*dt + exp(-X2*X2)*0.25*dW2;
")

mymod <- pomp(data=mymod.dat,time="t",t0=0,
              rprocess=euler.sim(step.fun=step.fun,delta.t=myh),
              statenames=c("X1","X2"),paramnames=c("theta1","theta2","theta3"))

rmeas <- Csnippet("
 Y1 = X1 + rnorm(0,1e-2);
 Y2 = X2 + rnorm(0,1e-2);
")

mymod <- pomp(mymod,rmeasure=rmeas,statenames=c("X1","X2"))

dmeas <- Csnippet("
  lik = dnorm(Y1,X1,1e-2,1) + dnorm(Y2,X2,1e-2,1);
")

mymod <- pomp(mymod,dmeasure=dmeas,statenames=c("X1","X2"))

myprior <- function(z)
{   
    return(dnorm(x = z, mean = 0, sd = 100, log = TRUE))
}

myposterior <- function(den, prior)
{       
    den[den <= 2.2e-16] = 2.2e-16
    loglik = sum(log(den))
    return(loglik + sum(prior))
}

burnin = 100
numsteps = 20000
totsteps = numsteps + burnin
numparticles = 10

artrack = numeric(length=(totsteps-1))
x = matrix(0, nrow = totsteps, ncol = 3)
x[1,] = c(0.1,0.1,0.1)

oldpf <- pfilter(mymod, Np = numparticles, params = c(theta1 = x[1,1], theta2 = x[1,2], theta3 = x[1,3], X1.0 = mymod.dat[1,]$Y1, X2.0 = mymod.dat[1,]$Y2))
oldden <- logLik(oldpf)
oldpost <- myposterior(oldden, myprior(x[1,]))

for(i in c(1:(totsteps-1))) {
        z = rnorm(n = 3, sd = 0.1)
        prop = x[i,] + z

        proppf <- pfilter(mymod, Np = numparticles, params = c(theta1 = prop[1], theta2 = prop[2], theta3 = prop[3], X1.0 = mymod.dat[1,]$Y1, X2.0 = mymod.dat[1,]$Y2))
        propden <- logLik(proppf)
        proppost <- myposterior(propden, myprior(prop))
        # print(propden)

        rho = exp(propden - oldden)

        u = runif(n = 1)
        if(rho > u) {
                x[i+1,] = prop
                oldpf = proppf
                oldden = propden
                oldpost = proppost
                artrack[i+1] = 1
                # print(paste("Accepted step", i, ": ", paste("theta[", c(1:2), "]=", format(prop, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        }
        else {
                x[i+1,] = x[i,]
                artrack[i+1] = 0
                # print(paste("Rejected step", i, ": ", paste("theta[", c(1:2), "]=", format(prop, digits = 3, scientific = TRUE), collapse = ', ', sep = '')))
        }
        if ((i %% 100) == 0) print(i)
        flush.console()
}

finaltime = proc.time() - ptm
print(finaltime)

artrack = artrack[(burnin+1):totsteps]
arratio = sum(artrack)/length(artrack)

save.image(file="pomp_vanderpol1.RData")
