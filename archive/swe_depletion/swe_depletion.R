##########################
#
# SWE depletion analyses 
# for picker project
#
# MML 04/24/18
##########################

rm(list=ls())

require(raster)
require(ncdf4)
require(rgdal)
require(gdalUtils)
require(lubridate)

setwd("L:/data_repo/gis_data/")

# define the projection - EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 

########## SWE data from Mudryk ##########
# swe data in m from Mudryk
# daily mean SWE (GS2,MERRA2,Brown,Crocus)
# list the monthly files
swe.files <- list.files(pattern =".nc",path ="swe_mudryk_blended/",full.names=T)

# read one file in to use for reprojecting
pr <- raster(swe.files[1])
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')

########## MODIS tree cover ##########
vcf <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')
vcf <- projectRaster(vcf,pr,overwrite=T,
                     filename = "L:/projects/boreal_swe_depletion/data/MOD44B_2014_mosaic_50km_ease.tif")

########## ETOPO1 Elevation ##########
#https://www.ngdc.noaa.gov/mgg/global/global.html
# info on grid vs. cell registration
topo <- raster("etopo1/original/ETOPO1_Ice_c_geotiff.tif")
projection(topo) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# aggregate to 0.5 degree resolution and calculate st dev of elevation
# use to identify mountains (Mudryk et al, 2017)
topo <- aggregate(topo,fact=30,fun=sd,overwrite=T,
                   filename="L:/projects/boreal_swe_depletion/data/ETOPO1_Ice_c_0.5deg_stdev.tif")
#reproject to EASE 2.0
topo <- projectRaster(topo,pr,overwrite=T,
                      filename = "L:/projects/boreal_swe_depletion/data/ETOPO1_Ice_c_50km_ease_stdev.tif")
# reclassify to create mountain mask
topo <- reclassify(topo,c(0,200,1,200,Inf,NA),overwrite=T,
                   filename="L:/projects/boreal_swe_depletion/data/ETOPO1_Ice_c_50km_ease_mtn_mask.tif")

########## GLC2000 Land Cover ##########
glc <- raster("GLC2000/original/glc2000_v1_1_Tiff/Tiff/glc2000_v1_1.tif")
projection(glc) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
glc <- crop(glc,c(-180,180,50,90))
# glc.laea <- projectRaster(glc,crs=laea,res=1,
#                           filename="L:/projects/boreal_swe_depletion/data/GLC2000_laea_1km.tif")

glc.legend <- read.csv('GLC2000/original/glc2000_v1_1_Tiff/Tiff/Global_Legend.csv',header=T)[,1:2]
## write function to calculate mode, for resampling
get.mode <- function(x,na.rm=T){
  f <- table(x)
  as.numeric(names(f)[which.max(f)])
}

get.mode2 <- function(x,na.rm=T){
  f <- table(x)
  as.numeric(names(sort(f,decreasing=T))[2])
}

## write function to calculate frequency of mode, for resampling
get.mode.freq <- function(x,na.rm=T){
  f <- table(x)
  (f)[which.max(f)]
}

get.mode2.freq <- function(x,na.rm=T){
  f <- table(x)
  sort(f,decreasing=T)[2]
}

## aggregate GLC2000 to ~0.5 degree resolution using mode, and make a map of the frequency
glc.mode <- aggregate(glc,fact=56,fun=get.mode,progress="text",overwrite=T,
                      filename="L:/projects/boreal_swe_depletion/data/GLC2000_0.5_mode.tif")  
glc.mode.freq <- aggregate(glc,fact=56,fun=get.mode.freq,progress="text",overwrite=T,
                           filename="L:/projects/boreal_swe_depletion/data/GLC2000_0.5_mode_freq.tif") 

glc.mode2 <- aggregate(glc,fact=56,fun=get.mode2,progress="text",overwrite=T,
                      filename="L:/projects/boreal_swe_depletion/data/GLC2000_0.5_mode2.tif")  
glc.mode2.freq <- aggregate(glc,fact=56,fun=get.mode2.freq,progress="text",overwrite=T,
                           filename="L:/projects/boreal_swe_depletion/data/GLC2000_0.5_mode_freq2.tif") 

## reproject aggregated GLC2000 to EASE grid resolution
glc.mode <- projectRaster(glc.mode,pr,method='ngb',progress='TEXT',overwrite=T,
                          filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode.tif')

glc.mode.freq <- projectRaster(glc.mode.freq,pr,method='ngb',progress='TEXT',overwrite=T,
                               filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode_freq.tif')

glc.mode2 <- projectRaster(glc.mode2,pr,method='ngb',progress='TEXT',overwrite=T,
                          filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode2.tif')

glc.mode2.freq <- projectRaster(glc.mode2.freq,pr,method='ngb',progress='TEXT',overwrite=T,
                               filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode2_freq.tif')

# ## create mask files for 50% and 75% vegetation coverage in each class
# glc.ease.50.mask <- reclassify(glc.mode.freq,matrix(c(0,1567,NA,1568,3136,50),ncol=3,byrow=T),right=NA,overwrite=T,
#                                filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode_freq_50.tif')
# 
# glc.ease.75.mask <- reclassify(glc.mode.freq,matrix(c(0,2351,NA,2352,3136,75),ncol=3,byrow=T),right=NA,overwrite=T,
#                                filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode_freq_75.tif')
# 
# ## mask glc to include majority pixels only
# glc.50pct.ease <- mask(glc.ease,glc.ease.50.mask,maskvalue=50,updatevalue=NA,inverse=T,overwrite=T,
#                        filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.25_50_mask.tif')
# 
# glc.75pct.ease <- mask(glc.ease,glc.ease.75.mask,maskvalue=75,updatevalue=NA,inverse=T,overwrite=T,
#                        filename='L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.25_50_mask.tif')




###################################
###################################
## load all of the raw data sets ##
###################################
###################################
topo <- raster("L:/projects/boreal_swe_depletion/data/ETOPO1_Ice_c_50km_ease_mtn_mask.tif")
glc.mode <- raster('L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode.tif')
glc.mode.freq <- raster('L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode_freq.tif')
glc.mode2 <- raster('L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode2.tif')
glc.mode2.freq <- raster('L:/projects/boreal_swe_depletion/data/GLC2000_EASE_from_0.5_mode2_freq.tif')
vcf <- raster("L:/projects/boreal_swe_depletion/data/MOD44B_2014_mosaic_50km_ease.tif")

# swe data in m from Mudryk
# daily mean SWE (GS2,MERRA2,Brown,Crocus)
# list the monthly files
swe.files <- list.files(pattern =".nc",path ="swe_mudryk_blended/",full.names=T)
era.files <- list.files(pattern = ".nc", path = "era_interim/air_temp_2m_spring/",full.names=T)

## read in file from first year
swe <- stack(swe.files[1])
## get date
ts <- swe@z
ts <- strptime(ts$Date,"%Y-%m-%d",tz="GMT")
d <- cbind(year(ts),month(ts),julian(ts,origin=ts[1])+1)
d <- d[which(d[,2]<7 & d[,2]>1),]
# change the year, just for 1981 - becuase there appears to be an error
d[,1] <- 1981
#subset to Feb-June
swe <- subset(swe,which(as.numeric(substr(names(swe),7,8))>=2 & as.numeric(substr(names(swe),7,8))<7))
# crop to boreal region
swe <- crop(swe,c(-180,180,50,90))
# reproject
swe <- projectRaster(swe,pr)
# reclass to set  values less than 1cm to NA
swe <- reclassify(swe,rcl=c(-Inf,0.01,NA))

# era interim temp data
era <- brick(era.files[1])
#reproject era
era <- projectRaster(era,pr)
#calculate daily means
era.d <- overlay(subset(era,seq(1,nlayers(era),4)),
             subset(era,seq(2,nlayers(era),4)),
             subset(era,seq(3,nlayers(era),4)),
             subset(era,seq(1,nlayers(era),4)),
             fun=mean)


# put all of the data together for the model
mod.dat <- cbind(rep(1:ncell(swe),nlayers(swe)),
                 rep(coordinates(swe)[,1],nlayers(swe)),
                 rep(coordinates(swe)[,2],nlayers(swe)),
                 rep(d[,1],each=ncell(swe)),
                 rep(d[,3],each=ncell(swe)),
                 as.vector(getValues(swe)),
                 as.vector(getValues(era.d)),
                 rep(getValues(vcf),nlayers(swe)),
                 rep(getValues(glc.mode),nlayers(swe)),
                 rep(getValues(glc.mode.freq),nlayers(swe)),
                 rep(getValues(glc.mode2),nlayers(swe)),
                 rep(getValues(glc.mode2.freq),nlayers(swe)))
#                 rep(getValues(topo),nlayers(swe)))
colnames(mod.dat) <- c("cell","x.coord","y.coord","year","jday","swe","t.air",
                        "vcf","glc1","glc1f","glc2","glc2f")

# remove NA values
mod.dat <- na.omit(mod.dat)

# write data to csv file
write.csv(mod.dat,file="L:/projects/boreal_swe_depletion/data/swe_depletion_model_data_vcf_no_topo.csv",
          row.names=F)

removeTmpFiles(h=0)
########################################
#### now loop through all 35 years #####
########################################
for(i in 2:length(swe.files))
{
  swe <- brick(swe.files[i])
  ## get date
  ts <- swe@z
  ts <- strptime(ts$Date,"%Y-%m-%d",tz="GMT")
  d <- cbind(year(ts),month(ts),julian(ts,origin=ts[1])+1)
  d <- d[which(d[,2]<7 & d[,2]>1),]

  #subset to Feb-June
  swe <- subset(swe,which(as.numeric(substr(names(swe),7,8))>=2 & as.numeric(substr(names(swe),7,8))<7))
  # crop to boreal region
  swe <- crop(swe,c(-180,180,50,90))
  # reproject
  swe <- projectRaster(swe,pr)
  # reclass to set  values less than 1cm to NA
  swe <- reclassify(swe,rcl=c(-Inf,0.01,NA))
  
  # era interim temp data
  era <- brick(era.files[i])
  #reproject era
  era <- projectRaster(era,pr)
  #calculate daily means
  era.d <- overlay(subset(era,seq(1,nlayers(era),4)),
                   subset(era,seq(2,nlayers(era),4)),
                   subset(era,seq(3,nlayers(era),4)),
                   subset(era,seq(1,nlayers(era),4)),
                   fun=mean)
  
  
  # put all of the data together for the model
  new.dat <- cbind(rep(1:ncell(swe),nlayers(swe)),
                   rep(coordinates(swe)[,1],nlayers(swe)),
                   rep(coordinates(swe)[,2],nlayers(swe)),
                   rep(d[,1],each=ncell(swe)),
                   rep(d[,3],each=ncell(swe)),
                   as.vector(getValues(swe)),
                   as.vector(getValues(era.d)),
                   rep(getValues(vcf),nlayers(swe)),
                   rep(getValues(glc.mode),nlayers(swe)),
                   rep(getValues(glc.mode.freq),nlayers(swe)),
                   rep(getValues(glc.mode2),nlayers(swe)),
                   rep(getValues(glc.mode2.freq),nlayers(swe)))
#                   rep(getValues(topo),nlayers(swe)))

  mod.dat <- rbind(mod.dat,new.dat)
  # remove NA values
  mod.dat <- na.omit(mod.dat)
  
  # write data to csv file
  write.csv(mod.dat,file="L:/projects/boreal_swe_depletion/data/swe_depletion_model_data_vcf_no_topo.csv",
            row.names=F,append=T)
  rm(swe,era,era.d)
  removeTmpFiles(h=0)
}

############## Archived Code ######################
###################################################
### load, reproject, and mosaic MODIS VCF data  ###
###################################################
vcf.2014.w <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_west_mosaic_sin.tif')
vcf.2014.e <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_east_mosaic_sin.tif')

r <- matrix(c(101,255,NA),ncol=3)

vcf.2014.w <- reclassify(vcf.2014.w,r,progress='text')

vcf.2014.e <- reclassify(vcf.2014.e,r,progress='text')

vcf.2014.w.ease <- projectRaster(vcf.2014.w,pr,progress='text')
vcf.2014.e.ease <- projectRaster(vcf.2014.e,pr,progress='text')

vcf.2014.ease <- mosaic(vcf.2014.w.ease,vcf.2014.e.ease,fun=mean,overwrite=T,  
                        filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')

vcf <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')



