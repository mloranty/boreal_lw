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
pairs(~ analysisDF$melt.mm.day + 
        analysisDF$log.melt + 
        analysisDF$lat + 
        analysisDF$meltTempC + 
        analysisDF$maxSwe.m)

#calculate log Abs melt rate

datalist <- #jags regression
  datalist <- list(Nobs= dim(cellSwe7)[1],
                   b0=cellSwe7$logAbsRate,
                   glcIDB=cellSwe7$gcID,
                   TempAB=cellSwe7$tair,
                   CanopyB=cellSwe7$vcf,
                   SweMax= cellSwe7$sweMax,
                   sweDay=cellSwe7$meltStart,
                   GCyearB=cellSwe7$gcyearID,
                   Ngcyear=dim(epsTable)[1],
                   ygcIDB=epsTable$gcID,
                   startb=startID,
                   endb=endID,
                   Nglc=dim(IDSglc)[1],
                   cellID=cellSwe7$cellID,
                   Ncell=nrow(cellDF),
                   TempMean=tempMean,
                   CanopyMean=CanopyMean,
                   SdayMean=SdayMean,
                   MaxMean=MaxMean)
