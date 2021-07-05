library(raster)
require(ncdf4)
require(rgdal)
require(gdalUtils)
require(lubridate)
library(tmap)


# define the projection - EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
#EPSG: 6931


###########################################
# read in data
# list the monthly files
swe.files <- list.files(pattern =".nc",path ="E:/Google Drive/GIS/swe_mudryk_blended",full.names=TRUE)

# read one file in to use for reprojecting
pr <- raster(swe.files[1])

#https://www.ngdc.noaa.gov/mgg/global/global.html
# info on grid vs. cell registration
topo <- raster("E:/Google Drive/GIS/boreal_swe_all_data/ETOPO1_Ice_c_geotiff.tif")

#modis tree cover
vcf.ease <- raster("E:/Google Drive/GIS/boreal_swe_all_data/archive/data/MOD44B_2014_mosaic_50km_ease.tif")

#glc 2000
glc <- raster("E:/Google Drive/GIS/boreal_swe_all_data/glc2000_v1_1.tif")

#read in swe files
sweAll <- list(stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2000.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2001.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2002.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2003.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2004.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2005.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2006.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2007.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2008.nc"),
 stack("E:/Google Drive/GIS/swe_mudryk_blended/SWE_obsMEAN4x.2009.nc"))



########## SWE data from Mudryk ##########
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')

########## ETOPO1 Elevation ##########

topo@crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
topoN <- crop(topo, extent(-180,180,45,90))


topo.ease <- projectRaster(topoN, crs=laea)

topo.sd <- aggregate(topo.ease,fact=50000/res(topo.ease),fun=sd)


topo.mask <- reclassify(topo.sd, matrix(c(0,200,1,
                                        200.00001,3000,0),byrow=TRUE))

#final mask
topo.maskR <- resample(topo.mask, pr, method="ngb")

#test mask
pr.m <- mask(pr, topo.maskR,maskvalue=0)
 
########### GLC -----

#define projection
projection(glc) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 '
#crop
glc <- crop(glc,c(-180,180,50,90))

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
glc.mode <- aggregate(glc,fact=56,fun=get.mode)  
glc.mode.freq <- aggregate(glc,fact=56,fun=get.mode.freq) 

glc.mode2 <- aggregate(glc,fact=56,fun=get.mode2)  
glc.mode2.freq <- aggregate(glc,fact=56,fun=get.mode2.freq) 


## reproject aggregated GLC2000 to EASE grid resolution
glc.mode.ease <- projectRaster(glc.mode,pr,method='ngb')

glc.mode.freq.ease <- projectRaster(glc.mode.freq,pr,method='ngb')

glc.mode2.ease <- projectRaster(glc.mode2,pr,method='ngb')

glc.mode2.freq.ease <- projectRaster(glc.mode2.freq,pr,method='ngb')


#subset most frequent landclass of interest
#
#4: Tree Cover, needle-leaved, evergreen
#5: Tree Cover, needle-leaved, deciduous
#6: Tree Cover, mixed leaf type
#12: Shrub Cover, closed-open, deciduous
#13: Herbaceous Cover, closed-open

#create function
glcSub <- function(x){
  ifelse(x == 4 | x == 5 | x == 6 | x == 12 | x == 13, x, NA)
  
}


glc.reclass <- calc(glc.mode.ease, glcSub)


tm_shape(glc.reclass)+
  tm_raster(style="cat",palette="Pastel1", 
            labels=c("4: Tree Cover, needle-leaved, evergreen",
              "5: Tree Cover, needle-leaved, deciduous",
              "6: Tree Cover, mixed leaf type",
              "12: Shrub Cover, closed-open, deciduous",
              "13: Herbaceous Cover, closed-open"))+
  tm_layout(legend.outside = TRUE)


#calculate proportion of most frequent glc
glc.mode.p.ease <- glc.mode.freq.ease/3136
glc.mode2.p.ease <- glc.mode2.freq.ease/3136

plot(glc.mode.p.ease)
plot(glc.mode2.p.ease)
#takes only majority land cover
glcP.mask <- reclassify(glc.mode.p.ease, matrix(c(0,0.5,NA,
                                                  0.5,1,1), byrow=TRUE, ncol=3))
plot(glcP.mask)

#get only majority land cover
glc.maj <- mask(glc.reclass,glcP.mask)
#get land cover not in mnts
glc.maj2 <- mask(glc.maj, topo.maskR)

tm_shape(glc.maj2)+
  tm_raster(style="cat",palette="Pastel1", colorNA = "grey90", 
            labels=c("4: Tree Cover, needle-leaved, evergreen",
                     "5: Tree Cover, needle-leaved, deciduous",
                     "6: Tree Cover, mixed leaf type",
                     "12: Shrub Cover, closed-open, deciduous",
                     "13: Herbaceous Cover, closed-open"))+
  tm_layout(legend.outside = TRUE, title = "Majority LC no mnts")

pr.m2 <- mask(pr.m, glcP.mask,maskvalue=NA)
plot(pr.m2)


######match all variables to masked area ----
vcf.mask <- mask(vcf,glc.maj2)


####### daily swe data ----
sweDates <- as.Date(names(sweAll[[1]]), "X%Y.%m.%d")
sweDOY <- yday(sweDates)
sweMonth <- month(sweDates)
plot(sweAll[[1]][[1]])
#subset stack to be only the melt period
swePeriod <-sweAll[[1]][[which(sweMonth >= 2 & sweMonth <=6)]]
plot(swePeriod)
#work with all swe data
#start with just 1
sweA.ease <- projectRaster(swePeriod, pr)
#apply mask for glc and land cover
sweA.mask <- mask(sweA.ease, glc.maj2)

##########################
##### Filter point  #####
##########################
#get max and min throughout the period

sweA.Max <- calc(sweA.mask, fun=max, na.rm=TRUE )
sweA.Min <- calc(sweA.mask, fun=min, na.rm=TRUE )


#exclude sites that don't get over 4 cm of swe
#swe units in meters
#adds to glc.maj2 since those are already excluded from swe
max.thresh <- function(x){
  ifelse(x <= 0.04, NA, x)
}

yearMask1 <- calc(sweA.Max, fun=max.thresh )

#now apply mask to Max, Min and sweA.mask
sweA.mask2 <- mask(sweA.mask, yearMask1)
sweMax.mask <- mask(sweA.Max, yearMask1)
sweMin.mask <- mask(sweA.Min, yearMask1)

#get last day of within 80 % of max
thrsh80 <- function(x, y){
  ifelse(x >= 0.8*y,1,NA )
}
#get function if swe is within the 80% threshold
swe80 <- overlay(sweA.mask2,sweMax.mask, fun=thrsh80)
#put in table
swe80Table <- getValues(swe80)
#get doy for swe in data set
sweADoy <- yday(as.Date(names(swePeriod), "X%Y.%m.%d"))
#function to get last day in column equal to 1
layer80 <- function(x){
 max(which(x == 1))
}
#apply to table
swe80Layer <- apply(swe80Table,1,layer80)
swe80LayerA <- ifelse(swe80Layer == -Inf,NA,swe80Layer)
#add back into the raster
meltStartDay <- setValues(sweMax.mask,swe80LayerA )
plot(meltStartDay)
#get actual swe value on start of melt day


##??? need this? also currently don't need original filter 4
# reclass to set  values less than 1cm to NA
swe <- reclassify(swe,rcl=c(-Inf,0.01,NA))


tm_shape(sweMax.mask)+
  tm_raster(title="2000 Maximum swe", palette= "BuPu",style="quantile")+
  tm_layout(legend.outside=TRUE)
#create a mask for snow extent


