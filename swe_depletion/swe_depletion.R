
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

setwd("L:/data_repo/gis_data/")

###################################
###################################
## load all of the raw data sets ##
###################################
###################################


#list the monthly files
swe.files <- list.files(pattern='.nc',path='swe_mudryk_blended/',full.names=T)

test <- stack(swe.files[1])
# read data for March, 1980-2014
# make a Raster stack from a list of file paths - subset the filenames to include only months of 

## the following is not necessary for the swe decline, but has more examples

################################################
###  load era temp data for March and April  ###
################################################
era <- brick("era_interim/original/interim_1981-02-01to1981-06-30_00061218.nc")
################################################
###  load CRU temp data for March and April  ###
################################################

cru.tmp <- brick('cru/cru_ts4.01/cru_ts4.01.1901.2016.tmp.dat.nc')
cru.tmp <- crop(cru.tmp,c(-180,180,50,90))

# reproject to match GlobSnow #
cru.ease <- projectRaster(cru.tmp,swe.dif.pos,filename='cru/cru_ts3.24.1901.2015.tmp.dat.ease.N.tif',overwrite=T)

# subset March data for 1980-2014 #
cru.mar <- subset(cru.ease,
                  which(as.numeric(substr(cru.tmp@z$time,1,4))>1979 &
                          as.numeric(substr(cru.tmp@z$time,1,4))<2015 &
                          as.numeric(substr(cru.tmp@z$time,6,7))==3))

# subset April data for 1980-2014 #
cru.apr <- subset(cru.ease,
                  which(as.numeric(substr(cru.tmp@z$time,1,4))>1979 &
                          as.numeric(substr(cru.tmp@z$time,1,4))<2015 &
                          as.numeric(substr(cru.tmp@z$time,6,7))==4))

# subset May data for 1980-2014 #
cru.may <- subset(cru.ease,
                  which(as.numeric(substr(cru.tmp@z$time,1,4))>1979 &
                          as.numeric(substr(cru.tmp@z$time,1,4))<2015 &
                          as.numeric(substr(cru.tmp@z$time,6,7))==5))

writeRaster(cru.mar,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Mar.ease.N.tif')
writeRaster(cru.apr,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Apr.ease.N.tif')
writeRaster(cru.may,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.May.ease.N.tif')

# calculate the difference #
cru.dif <- cru.apr-cru.mar
writeRaster(cru.dif,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Mar_Apr.dif.ease.N.tif')

# reclass to include only positive values
cru.dif.pos <- reclassify(cru.dif,r,overwrite=T,
                          filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Apr_Mar_pos.ease.N.tif')

# calculate swe reduction per unit temp increase using positive datasets
swe.cru <- swe.dif.pos/cru.dif.pos
writeRaster(swe.cru,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/swe_cru_Mar_Apr_change.tif')

###################################################
### load, reproject, and mosaic MODIS VCF data  ###
###################################################
vcf.2014.w <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_west_mosaic_sin.tif')
vcf.2014.e <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_east_mosaic_sin.tif')

r <- matrix(c(101,255,NA),ncol=3)

vcf.2014.w <- reclassify(vcf.2014.w,r,progress='text', 
                         filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_west_mosaic_sin_rcl.tif')

vcf.2014.w.ease <- projectRaster(vcf.2014.w,cru.dif,progress='text',
                                 filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_west_mosaic_ease.tif')

vcf.2014.e <- reclassify(vcf.2014.e,r,progress='text',overwrite=T,  
                         filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_east_mosaic_sin_rcl.tif')

vcf.2014.e.ease <- projectRaster(vcf.2014.e,cru.dif,progress='text',overwrite=T,  
                                 filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_east_mosaic_ease.tif')

vcf.2014.ease <- mosaic(vcf.2014.w.ease,vcf.2014.e.ease,fun=mean,overwrite=T,  
                        filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')

#vcf.2014.ease <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')

###################################################
### load daily SWE files from Lawrence Mudryk    ##
### these are a combo of GlobSnow, MERRA & Brown ##
###################################################

###SKIP FOR NOW - HAVING PROBLEMS WITH EASE 2.0 GRID PROJECTION###
#setwd("Y:/swe_mudryk/")
swe <- brick('SWE_obsMEAN_M2BGS_ease2.2010.nc')

## note this is the correct PROJ.4 description for EASE2, but not sure this works with Lawrence's files ##
projection(swe) <- "+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m "


###################################################
### load GLC2000 landcover map                   ##
###                                              ##
###################################################
glc <- raster('GLC2000/Tiff/glc2000_v1_1.tif')
projection(glc) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
glc <- crop(glc,c(-180,180,50,90))
glc.legend <- read.xlsx('/Users/mloranty/Google Drive/GIS_Data/GLC2000/Tiff/Global_Legend.xls',sheetIndex=1)[,1:2]
colnames(glc.legend) <- c('zone','names')

#pr <- projectRaster(vcf.2014.ease,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",method='ngb')

## write function to calculate mode, for resampling
get.mode <- function(x,na.rm=T){
  f <- table(x)
  as.numeric(names(f)[which.max(f)])
}

## write function to calculate frequency of mode, for resampling
get.mode.freq <- function(x,na.rm=T){
  f <- table(x)
  (f)[which.max(f)]
}

## aggregate GLC2000 to ~0.25 degree resolution using mode, and make a map of the frequency
glc.mode <- aggregate(glc,fact=25,fun=get.mode)  
glc.mode.freq <- aggregate(glc,fact=25,fun=get.mode.freq) 

## reproject aggregated GLC2000 to EASE grid resolution
glc.ease <- projectRaster(glc.mode,cru.ease,method='ngb',progress='TEXT',overwrite=T,
                          filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode.tif')

glc.ease.freq <- projectRaster(glc.mode.freq,cru.ease,method='ngb',progress='TEXT',overwrite=T,
                               filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode_freq.tif')

## create mask files for 50% and 75% vegetation coverage in each class
glc.ease.50.mask <- reclassify(glc.ease.freq,matrix(c(0,312,NA,313,625,50),ncol=3,byrow=T),right=NA,overwrite=T,
                               filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode_freq_50.tif')

glc.ease.75.mask <- reclassify(glc.ease.freq,matrix(c(0,468,NA,469,615,75),ncol=3,byrow=T),right=NA,overwrite=T,
                               filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode_freq_75.tif')

## mask glc to include majority pixels only
glc.50pct.ease <- mask(glc.ease,glc.ease.50.mask,maskvalue=50,updatevalue=NA,inverse=T,overwrite=T,
                       filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_50_mask.tif')

glc.75pct.ease <- mask(glc.ease,glc.ease.75.mask,maskvalue=75,updatevalue=NA,inverse=T,overwrite=T,
                       filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_50_mask.tif')

## STOPPING HERE FOR NOW ##
dat <- ls()
save.image(dat,file="C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/monthly_swe_analysis_23Nov.RData")

#############################################
#############################################
##                                         ##
## LOAD ALL PREVIOUSLY PROCESSED DATA SETS ##
##                                         ##
#############################################
#############################################



