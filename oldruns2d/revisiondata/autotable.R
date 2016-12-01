rm(list=ls(all=TRUE))

library('stargazer')

options(digits = 3)
truemean = c(1, 1, 4)

rmse <- function(app, tru)
{
    z = app/tru - 1
    return(mean(abs(z)))
}

compstats <- function(ms, ar)
{
    meansamp = apply(ms^2, 2, mean)
    sdsamp = apply(ms^2, 2, sd)
    arsamp = mean(ar)
    errsamp = rmse(meansamp, truemean)
    return(c(meansamp=meansamp,sdsamp=sdsamp,arsamp=arsamp,errsamp=errsamp))
}

fnames = c('samples_pomp_noadapt_vanderpol_by2.RData',
'samples_pomp_noadapt_vanderpol_by4.RData',
'samples_pomp_noadapt_vanderpol_by8.RData',
'samples_pomp_vanderpol_by2.RData',
'samples_pomp_vanderpol_by4.RData',
'samples_pomp_vanderpol_by8.RData',
'samples_vanderpol_by1.RData',
'samples_vanderpol_by2.RData',
'samples_vanderpol_by4.RData',
'samples_vanderpol_by8.RData')

mymat = matrix(nrow=length(fnames),ncol=9)

for (iii in c(1:length(fnames)))
{
    load(fnames[iii])
    mymat[iii,9] = myh
    mymat[iii,1:8] = compstats(mcmcsamples[1:10000,], artrack)
}

type = c(rep('pomp',3),rep('pomp-adaptive',3),'Eulerian',rep('DTQ',3))

mydf = data.frame(mymat,type)
names(mydf)[1:9] = c('theta1','theta2','theta3','sd1','sd2','sd3','ar','err','h')

# first table that includes only non-adaptive MCMC results
stargazer(mydf[c(c(7:10),c(1:3)),-c(4:6)],summary=FALSE)

# second table with adaptive MCMC results for pomp only
stargazer(mydf[c(4:6),-c(4:6)],summary=FALSE)






