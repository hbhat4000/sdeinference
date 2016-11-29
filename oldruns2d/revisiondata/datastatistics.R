rm(list=ls(all=TRUE))

library('Metrics')

options(digits = 3)
truemean = c(1, 1, 4)

load('samples_pomp_vanderpol_by2.RData')
pomp2 = mcmcsamples
arpomp2 = mean(artrack)
meanpomp2 = c(mean(pomp2[,1]^2), mean(pomp2[,2]^2), mean(pomp2[,3]^2))
diffmeanpomp2 = truemean - meanpomp2
sdpomp2 = c(sd(pomp2[,1]^2), sd(pomp2[,2]^2), sd(pomp2[,3]^2))
rmsepomp2 = rmse(meanpomp2, truemean)

load('samples_pomp_vanderpol_by4.RData')
pomp4 = mcmcsamples
arpomp4 = mean(artrack)
meanpomp4 = c(mean(pomp4[,1]^2), mean(pomp4[,2]^2), mean(pomp4[,3]^2))
diffmeanpomp4 = truemean - meanpomp4
sdpomp4 = c(sd(pomp4[,1]^2), sd(pomp4[,2]^2), sd(pomp4[,3]^2))
rmsepomp4 = rmse(meanpomp4, truemean)

load('samples_pomp_vanderpol_by8.RData')
pomp8 = mcmcsamples
arpomp8 = mean(artrack)
meanpomp8 = c(mean(pomp8[,1]^2), mean(pomp8[,2]^2), mean(pomp8[,3]^2))
diffmeanpomp8 = truemean - meanpomp8
sdpomp8 = c(sd(pomp8[,1]^2), sd(pomp8[,2]^2), sd(pomp8[,3]^2))
rmsepomp8 = rmse(meanpomp8, truemean)

load('samples_vanderpol_by1.RData')
dtq1 = mcmcsamples
ardtq1 = mean(artrack)
meandtq1 = c(mean(dtq1[,1]^2), mean(dtq1[,2]^2), mean(dtq1[,3]^2))
diffmeandtq1 = truemean - meandtq1
sddtq1 = c(sd(dtq1[,1]^2), sd(dtq1[,2]^2), sd(dtq1[,3]^2))
rmsedtq1 = rmse(meandtq1, truemean)

load('samples_vanderpol_by2.RData')
dtq2 = mcmcsamples
ardtq2 = mean(artrack)
meandtq2 = c(mean(dtq2[,1]^2), mean(dtq2[,2]^2), mean(dtq2[,3]^2))
diffmeandtq2 = truemean - meandtq2
sddtq2 = c(sd(dtq2[,1]^2), sd(dtq2[,2]^2), sd(dtq2[,3]^2))
rmsedtq2 = rmse(meandtq2, truemean)

load('samples_vanderpol_by4.RData')
dtq4 = mcmcsamples
ardtq4 = mean(artrack)
meandtq4 = c(mean(dtq4[,1]^2), mean(dtq4[,2]^2), mean(dtq4[,3]^2))
diffmeandtq4 = truemean - meandtq4
sddtq4 = c(sd(dtq4[,1]^2), sd(dtq4[,2]^2), sd(dtq4[,3]^2))
rmsedtq4 = rmse(meandtq4, truemean)

load('samples_vanderpol_by8.RData')
dtq8 = mcmcsamples
ardtq8 = mean(artrack)
meandtq8 = c(mean(dtq8[,1]^2), mean(dtq8[,2]^2), mean(dtq8[,3]^2))
diffmeandtq8 = truemean - meandtq8
sddtq8 = c(sd(dtq8[,1]^2), sd(dtq8[,2]^2), sd(dtq8[,3]^2))
rmsedtq8 = rmse(meandtq8, truemean)

rmsevec = c(rmsepomp2, rmsepomp4, rmsepomp8, rmsedtq1, rmsedtq2, rmsedtq4, rmsedtq8)
artrackvec = c(arpomp2, arpomp4, arpomp8, ardtq1, ardtq2, ardtq4, ardtq8)