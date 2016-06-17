#clear memory
rm(list=ls(all=TRUE))
library('ggplot2')

#pick which directories and which variables we want
mydirs=c(6:7)
myvbls=c(1,2,255)

#colors
mycols=c("black","blue")
mypars=c("h=0.02","h=0.01")

#names of parameters
library(latex2exp)
cnames=c("theta1","theta2","sigeps2")
pnames=c(TeX('$\\theta_1$'),TeX('$\\theta_2$'),TeX('$\\log_{10}(\\sigma_{\\epsilon^2})$'))

#true values of parameters
tvals=c(0.5,1,log(0.01)/log(10))

#set burn-in period
burnin=c(1:100)

#log transform?
logtrans=c(FALSE,FALSE,TRUE)

#load data
mydata=list(NULL)
basedir="actual_bigrun"
for (j in c(1:length(mydirs)))
{
    jj=mydirs[j]
    fname=paste(basedir,jj,"/mcmc.out",sep='')
    temp=read.csv(file=fname,header=FALSE)
    mydata[[j]]=temp[,myvbls]
    mydata[[j]]=mydata[[j]][-burnin,]
    colnames(mydata[[j]])=cnames

    #apply log transform as needed
    for (k in c(1:length(myvbls)))
        if (logtrans[k])
            mydata[[j]][,k] = log(mydata[[j]][,k])/log(10)
}

#create intermediate data frames, more suitable for ggplot2
vbf=list(NULL)
for (k in c(1:length(myvbls)))
{
    for (j in c(1:length(mydirs)))
    {
        temp=data.frame(mydata[[j]][,k],mypar=rep(mypars[j],nrow(mydata[[j]])),stringsAsFactors=FALSE)
        colnames(temp)[1]=cnames[k]
        if (j==1) vbf[[k]]=temp
        else vbf[[k]]=rbind(vbf[[k]],temp)
    }
}

#Hadley's function
library('grid')
library('gridExtra')
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl$nrow=1
  grid.arrange(do.call(arrangeGrob,gl),legend,nrow = 2,heights = unit.c(unit(1, "npc") - lheight, lheight))
}

#create plots
myplots=list(NULL)
for (k in c(1:length(myvbls)))
{
    g=ggplot(vbf[[k]],aes_string(cnames[k]))
    g=g+stat_density(geom="line",aes(colour=mypar),bw="nrd",position="identity")
    g=g+theme_bw()
    g=g+xlab(pnames[k])
    g=g+geom_vline(xintercept=tvals[k],colour="red")
    g=g+scale_colour_manual("DTQ step",values=mycols,breaks=mypars)
    myplots[[k]]=g
}

#save plot
g=grid_arrange_shared_legend(myplots[[1]],myplots[[2]],myplots[[3]])
ggsave(filename="post_nonequi.pdf",plot=g,height=3,width=6,units="in")
ggsave(filename="post_nonequi.eps",plot=g,height=3,width=6,units="in")


