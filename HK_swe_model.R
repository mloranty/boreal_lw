###########################################################
###########################################################
###########################################################
############ SWE stats model v2.0              ############  
############ modified from v1.0                ############
############ H. Kropp                          ############
###########################################################
###########################################################
###########################################################


###########################################
########## Get melt data        ----- 

#read in swe data
source("c:/Users/hkropp/Documents/GitHub/boreal_lw/HK_swe_organize.R")

###########################################
########## Libraries        ----- 
library(raster)
library(ncdf4)
library(rgdal)
library(gdalUtils)
library(sp)
library(rjags)
library(coda)
library(mcmcplots)
library(loo)
library(dplyr)
###########################################
########## Directories         ----- 

#model out directory

modDir <- "E:/Google Drive/research/projects/boreal_swe/boreal_2021/model/run1"


###########################################
########## Set up model data        -----

#calculate Abs melt rate, all values represent decrease in swe
analysisDF$abs.melt <- abs(analysisDF$melt.mm.day)

#log transform
analysisDF$log.melt <- log(analysisDF$abs.melt)
#check covariate correlation
pairs(formula= ~ melt.mm.day + 
        log.melt + 
        vcf +
        lat + 
        meltTempC + 
        maxSwe.m, data=analysisDF)

pairs(formula= ~ melt.mm.day + 
        log.melt + 
        vcf +
        lat + 
        meltTempC + 
        maxSwe.m, data=analysisDF[glcID == 4,])

pairs(formula= ~ melt.mm.day + 
        log.melt + 
        vcf +
        lat + 
        meltTempC + 
        maxSwe.m, data=analysisDF[glcID == 5,])

pairs(formula= ~ melt.mm.day + 
        log.melt + 
        vcf +
        lat + 
        meltTempC + 
        maxSwe.m, data=analysisDF[glcID == 6,])

pairs(formula= ~ melt.mm.day + 
        log.melt + 
        vcf +
        lat + 
        meltTempC + 
        maxSwe.m, data=analysisDF[glcID == 12,])

pairs(formula= ~ melt.mm.day + 
        log.melt + 
        vcf +
        lat + 
        meltTempC + 
        maxSwe.m, data=analysisDF[glcID == 13,])

#check correlation
cor(analysisDF)

#create gcID column in table
colnames(glcID) <- c("glc","Desc")
glcID$gcID <- seq(1, nrow(glcID))

#join into analysis DF
analysisDFm1 <- left_join(analysisDF, glcID, by="glc")


#random effects ids
#need to organize table for eps ids
epsTable <- unique(data.frame(gcID=analysisDFm1$gcID,year=analysisDFm1$year))
epsTable <- epsTable[order(epsTable$gcID,epsTable$year),]
#this order will be by GCID
epsTable$gcyearID <- seq(1,dim(epsTable)[1])

#create index for averaging eps
gcIndT <- unique(data.frame(gcID=epsTable$gcID))
startID <- numeric(0)
endID <- numeric(0)

for(i in 1:dim(gcIndT)[1]){
  startID[i] <- head(which(epsTable$gcID==gcIndT$gcID[i]))[1]
  endID [i] <- tail(which(epsTable$gcID==gcIndT$gcID[i]))[6]
}

#join back into analysis DF
analysisDFm1 <- left_join(analysisDFm1, epsTable, by=c("gcID","year"))

#set up data for plotting
tempMean <- seq(floor(range(analysisDFm1$meltTempC)[1]),ceiling(range(analysisDFm1$meltTempC)[2]), length.out=200)
CanopyMean <- seq(floor(range(analysisDFm1$vcf)[1]),ceiling(range(analysisDFm1$vcf)[2]), length.out=200)
SdayMean <- seq(floor(range(analysisDFm1$doyStart)[1]),ceiling(range(analysisDFm1$doyStart)[2]), length.out=200)
MaxMean <- seq(0,0.5, length.out=200)




###########################################
########## Model        -----
#jags mixed effects regression
#organize into list
  datalist <- list(Nobs= dim(analysisDFm1)[1],
                   b0=analysisDFm1$logAbsRate,
                   glcIDB=analysisDFm1$gcID,
                   TempAB=analysisDFm1$meltTempC,#centered around 0
                   CanopyB=analysisDFm1$vcf,#centered at 20
                   SweMax= analysisDFm1$maxSwe.m, #centered at 0.15
                   sweDay=analysisDFm1$doyStart, #centered at 107
                   #GCyearB=analysisDFm1$gcyearID,
                   #Ngcyear=dim(epsTable)[1],
                   Nglc=dim(glcID)[1],
                   #ygcIDB=epsTable$gcID,
                   #startb=startID,
                   #endb=endID,
                   TempMean=tempMean,
                   CanopyMean=CanopyMean,
                   SdayMean=SdayMean,
                   MaxMean=MaxMean)

  inits <- list(list(tau.eb=rep(1,dim(gcIndT)[1])),
                list(tau.eb=rep(1.4,dim(gcIndT)[1])),
                list(tau.eb=rep(2,dim(gcIndT)[1])))
  
  
  parms <- c("betaB0S","betaB1","betaB2","betaB3","betaB4",
             "mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3","mu.betaB4",
             "sig.B0","sig.B1","sig.B2","sig.B3","sig.B4","trB0",
             "rep.b0","Dsum","loglike",
             "mu.Temp","mu.Canopy","mu.Onset","mu.Max")
  
  curve.mod <- jags.model(file="c:\\Users\\hkropp\\Documents\\GitHub\\boreal_lw\\HK_model_code.r",
                          data=datalist,
                          n.adapt=10000,
                          n.chains=3)
  
  curve.sample <- coda.samples(curve.mod,variable.names=parms,n.iter=3000,thin=1)						
  
  mcmcplot(curve.sample, parms=c(
    "betaB0S","betaB1","betaB2","betaB3","betaB4",
    "mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3","mu.betaB4",
    "sig.B0","sig.B1","sig.B2","sig.B3","sig.B4",
    "sig.vB","eps.bS","tau.eb","sig.es"),dir=paste0(modDir,"\\history"))
  
  #model output							   
  mod.out <- summary(curve.sample)
  
  
  write.table(mod.out$statistics,paste0(modDir,"\\curve_mod_stats.csv"),
              sep=",",row.names=TRUE)
  write.table(mod.out$quantiles,paste0(modDir,"\\curve_mod_quant.csv"),
              sep=",",row.names=TRUE)	
  
  
  #coda output
  chain1<-as.matrix(curve.sample[[1]])
  write.table(chain1,paste0(modDir,"\\chain1_coda.csv"), sep=",")
  chain2<-as.matrix(curve.sample[[2]])
  write.table(chain2,paste0(modDir,"\\chain2_coda.csv"), sep=",")
  chain3<-as.matrix(curve.sample[[3]])
  write.table(chain3,paste0(modDir,"\\chain3_coda.csv"), sep=",")		 