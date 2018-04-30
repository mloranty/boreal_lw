
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
                     filename = "MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_50km_ease.tif")

########## ETOPO1 Elevation ##########
#https://www.ngdc.noaa.gov/mgg/global/global.html
# info on grid vs. cell registration
topo <- raster("etopo1/original/ETOPO1_Ice_c_geotiff.tif")
projection(topo) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# aggregate to 0.5 degree resolution and calculate st dev of elevation
# use to identify mountains (Mudryk et al, 2017)
topo <- aggregate(topo,fact=30,fun=sd,
                   filename="etopo1/modified/ETOPO1_Ice_c_0.5deg_stdev.tif")
#reproject to EASE 2.0
topo <- projectRaster(topo,pr,
                      filename = "etopo1/modified/ETOPO1_Ice_c_50km_ease_stdev.tif")
# reclassify to create mountain mask
topo <- reclassify(topo,c(0,200,NA,200,Inf,1),overwrite=T,
                   filename="etopo1/modified/ETOPO1_Ice_c_50km_ease_mtn_mask.tif")

########## GLC2000 Land Cover ##########
glc <- raster("GLC2000/original/glc2000_v1_1_Tiff/Tiff/glc2000_v1_1.tif")
projection(glc) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
glc <- crop(glc,c(-180,180,50,90))
#glc.legend <- read.xlsx('/GLC2000/original/glc2000_v1_1_Tiff/Tiff//Global_Legend.xls',sheetIndex=1)[,1:2]
#colnames(glc.legend) <- c('zone','names')

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

## aggregate GLC2000 to ~0.5 degree resolution using mode, and make a map of the frequency
glc.mode <- aggregate(glc,fact=56,fun=get.mode,progress="text",
                      filename="GLC2000/modified/GLC2000_0.5_mode.tif")  
glc.mode.freq <- aggregate(glc,fact=56,fun=get.mode.freq,progress="text",
                           filename="GLC2000/modified/GLC2000_0.5_mode.tif") 

## reproject aggregated GLC2000 to EASE grid resolution
glc.ease <- projectRaster(glc.mode,cru.ease,method='ngb',progress='TEXT',overwrite=T,
                          filename='GLC2000/modified/GLC2000_EASE_from_0.5_mode.tif')

glc.ease.freq <- projectRaster(glc.mode.freq,cru.ease,method='ngb',progress='TEXT',overwrite=T,
                               filename='GLC2000/modified/GLC2000_EASE_from_0.5_mode_freq.tif')

## create mask files for 50% and 75% vegetation coverage in each class
glc.ease.50.mask <- reclassify(glc.ease.freq,matrix(c(0,312,NA,313,625,50),ncol=3,byrow=T),right=NA,overwrite=T,
                               filename='GLC2000/modified/GLC2000_EASE_from_0.5_mode_freq_50.tif')

glc.ease.75.mask <- reclassify(glc.ease.freq,matrix(c(0,468,NA,469,615,75),ncol=3,byrow=T),right=NA,overwrite=T,
                               filename='GLC2000/modified/GLC2000_EASE_from_0.5_mode_freq_75.tif')

## mask glc to include majority pixels only
glc.50pct.ease <- mask(glc.ease,glc.ease.50.mask,maskvalue=50,updatevalue=NA,inverse=T,overwrite=T,
                       filename='GLC2000/modified/GLC2000_EASE_from_0.25_50_mask.tif')

glc.75pct.ease <- mask(glc.ease,glc.ease.75.mask,maskvalue=75,updatevalue=NA,inverse=T,overwrite=T,
                       filename='GLC2000/modified/GLC2000_EASE_from_0.25_50_mask.tif')




###################################
###################################
## load all of the raw data sets ##
###################################
###################################

# swe data in m from Mudryk
# daily mean SWE (GS2,MERRA2,Brown,Crocus)
# list the monthly files
swe.files <- list.files(pattern =".nc",path ="swe_mudryk_blended/",full.names=T)
era.files <- list/files(pattern = ".nc", path = "era_interim/original/",full.names=T)




t2 <- projectRaster(test[[1]],res=50000,crs=laea,progress='text')
# read data for March, 1980-2014
# make a Raster stack from a list of file paths - subset the filenames to include only months of 

## the following is not necessary for the swe decline, but has more examples

################################################
###  load era temp data for March and April  ###
################################################
era <- brick("era_interim/original/interim_1981-02-01to1981-06-30_00061218.nc")







###################################################
### load GLC2000 landcover map                   ##
###                                              ##
###################################################

############## Archived Code ######################
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

vcf <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')



