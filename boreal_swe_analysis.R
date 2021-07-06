############################
# 
# updated code to streamline
# analysis workflow and fix
# maksing issue
#
# 07/01/21
############################

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

########## read in data ------
# list subset of monthly SWE files that include the 2000-2009 study period
swe.files <- list.files(pattern =".nc",path ="swe_mudryk_blended/",full.names=T)[20:29]
# read one file in to use for reprojecting
pr <- raster(swe.files[1])

# list subset of era files that include the 2000-2009 study period
era.files <- list.files(pattern = ".nc", path = "era_interim/air_temp_2m_spring/",full.names=T)[20:29]

# read reprojected modis vcf data
vcf <- raster("L:/projects/boreal_swe_depletion/data/MOD44B_2014_mosaic_50km_ease.tif")

# topography data
topo <- raster("etopo1/original/ETOPO1_Ice_c_geotiff.tif")

#GLC2000 Land Cover Data
glc <- raster("GLC2000/original/glc2000_v1_1_Tiff/Tiff/glc2000_v1_1.tif")

########## create SWE template to use for resampling ------
# using SWE layer read in above
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')

########## preprocess ETOPO1 Elevation Mask ------
#https://www.ngdc.noaa.gov/mgg/global/global.html
# info on grid vs. cell registration

projection(topo) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# crop to boreal
topo <- crop(topo,c(-180,180,45,90))

#reproject to EASE 2.0
topo.ease <- projectRaster(topo,crs=laea)

# aggregate to 0.5 degree resolution and calculate st dev of elevation
# use to identify mountains (Mudryk et al, 2017)
topo.sd <- aggregate(topo.ease,fact=50000/res(topo.ease),fun=sd)

# reclassify to create mountain mask
topo.mask <- reclassify(topo.sd,matrix(c(0,200,1,200,Inf,NA),ncol = 3, byrow = TRUE))

#reproject mask to match swe data
topo.maskR <- resample(topo.mask,pr, method = "ngb")

# test mask swe by topo to remove mountains
#swe.mask <- mask(pr,topo.maskR)


########## preprocess GLC2000 Land Cover ------
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
glc.mode <- aggregate(glc,fact=56,fun=get.mode,progress="text")  
glc.mode.freq <- aggregate(glc,fact=56,fun=get.mode.freq,progress="text") 

glc.mode2 <- aggregate(glc,fact=56,fun=get.mode2,progress="text")  
glc.mode2.freq <- aggregate(glc,fact=56,fun=get.mode2.freq,progress="text") 

## reproject aggregated GLC2000 to EASE grid resolution
glc.mode.ease <- projectRaster(glc.mode,pr,method='ngb',progress="text")

glc.mode.freq.ease <- projectRaster(glc.mode.freq,pr,method='ngb',progress="text")

glc.mode2.ease <- projectRaster(glc.mode2,pr,method='ngb',progress="text")

glc.mode2.freq.ease <- projectRaster(glc.mode2.freq,pr,method='ngb',progress="text")

########## preprocess SWE data -----

## subset, resample, and filter swe data
swe <- stack(swe.files[1])

# subset to Feb-June using layer names
swe <- subset(swe,which(as.numeric(substr(names(swe),7,8))>=2 & as.numeric(substr(names(swe),7,8))<7))
# crop to boreal region
swe <- crop(swe,c(-180,180,50,90))
# reproject
swe <- projectRaster(swe,pr)
# reclass to set  values less than 1cm to NA
swe <- reclassify(swe,rcl=c(-Inf,0.01,NA))

########## preprocess ERA Interim data -----
for(i in 1:length(era.files))
{
  # era interim temp data
  era <- brick(era.files[i])

  # output filename
  f <- gsub("00061218.nc", "daily.grd", era.files[i])
  
  # calculate daily means
  era.d <- overlay(subset(era,seq(1,nlayers(era),4)),
                   subset(era,seq(2,nlayers(era),4)),
                   subset(era,seq(3,nlayers(era),4)),
                   subset(era,seq(1,nlayers(era),4)),
                   fun=mean, progress = "text")
  
  # assign names based on the date from original file
  names(era.d) <- unique(substr(names(era),2,11))
  
  # round the extent 
  # note that this eliminates warnings during reprojection and a swath of missing data near 180 deg lon
  extent(era.d) <- round(extent(era.d))
  
  # reproject and write to file
  era.d <- projectRaster(era.d,pr, progress = "text",
                         filename = f, overwrite = T)

}


