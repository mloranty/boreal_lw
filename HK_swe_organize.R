library(raster)
require(ncdf4)
require(rgdal)
require(gdalUtils)
require(lubridate)
library(tmap)
library(sp)

###########################################
########## define projection ------
#  EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
#EPSG: 6931



###########################################
########## read in data          ----- 


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

#legnth of stack
NYears <- 10

###########################################
########## SWE prep -------
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')


###########################################
########## ETOPO1 Elevation -------

# topo needs projection defined
topo@crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# crop to general area
topoN <- crop(topo, extent(-180,180,45,90))
# reproject
topo.ease <- projectRaster(topoN, crs=laea)
# aggregate and get sd
topo.sd <- aggregate(topo.ease,fact=50000/res(topo.ease),fun=sd)

#anything below 200 can be used 1= use
topo.mask <- reclassify(topo.sd, matrix(c(0,200,1,
                                        200.00001,3000,0),byrow=TRUE))

#resample to match swe data
topo.maskR <- resample(topo.mask, pr, method="ngb")

#mask any area with sd over 2000
pr.m <- mask(pr, topo.maskR,maskvalue=0)
 

###########################################
########## GLC -----

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



###########################################
########## VCF ----
vcf.mask <- mask(vcf,glc.maj2)



###########################################
########## SWE melt calculations ----
sweDates <- list()
sweDOY <- list()
sweMonth <- list()
swePeriod <- list()
sweA.ease <- list()
sweA.mask <- list()
sweA.mask2 <- list()
sweA.Max <- list()
sweA.Min<- list()
for(i in 1:NYears){
  #get date information from layer names
  sweDates[[i]] <- as.Date(names(sweAll[[i]]), "X%Y.%m.%d")
  sweDOY[[i]] <- yday(sweDates[[i]])
  sweMonth[[i]] <- month(sweDates[[i]])
  

  #subset stack to be only the melt period
  swePeriod[[i]] <-sweAll[[i]][[which(sweMonth[[i]] >= 2 & sweMonth[[i]] <=6)]]
  
  #work with all swe data
  #start with just 1
  sweA.ease[[i]] <- projectRaster(swePeriod[[i]], pr)
  #apply mask for glc and land cover
  sweA.mask[[i]] <- mask(sweA.ease[[i]], glc.maj2)
  
  # reclass to set  values less than 1cm to NA
  #since these values
  sweA.mask2[[i]] <- reclassify(sweA.mask[[i]],rcl=c(-Inf,0.01,NA))
  
  #get max and min throughout the period
  
  sweA.Max[[i]] <- calc(sweA.mask2[[i]], fun=max, na.rm=TRUE )
  sweA.Min[[i]] <- calc(sweA.mask2[[i]], fun=min, na.rm=TRUE )
}


##### Filter point 1:
##### SWE max must be greater than 0.04m in a cell

#set max values below or equal to 0.04 to an NA
max.thresh <- function(x){
  ifelse(x <= 0.04, NA, x)
}
#apply to raster and use as a mask
yearMask1 <- list()
sweA.mask3  <- list()
sweMax.mask  <- list()
sweMin.mask  <- list()
for(i in 1:NYears){
  #adds to glc.maj2 since those are already excluded from swe
  yearMask1[[i]] <- calc(sweA.Max[[i]], fun=max.thresh )
  
  #now apply mask to Max, Min and sweA.mask
  sweA.mask3[[i]] <- mask(sweA.mask2[[i]], yearMask1[[i]])
  sweMax.mask[[i]] <- mask(sweA.Max[[i]], yearMask1[[i]])
  sweMin.mask[[i]] <- mask(sweA.Min[[i]], yearMask1[[i]])
}


##### Filter point 2:
##### any sites less 30 days of swe observations (due to <0.01) 


#get length of observations without NA
NA.func <- function(x){
 length(na.omit(x)) 
}

obsTable <- list()
obsCount <- list()
obsRaster <- list()
obsMask <- list()
sweA.mask4 <- list()
sweMax.mask2 <- list()
sweMin.mask2 <- list()
for(i in 1:NYears){
  #get actual observation count
  obsTable[[i]] <- getValues(sweA.mask3[[i]])
  
  obsCount[[i]] <- apply(obsTable[[i]], 1, NA.func)

  #add count of number of observations
  obsRaster[[i]] <- setValues(sweMax.mask[[i]], obsCount[[i]])

  #create mask
  obsMask[[i]] <- reclassify(obsRaster[[i]], c(-1,29,NA))



  #mask 
  #now apply mask to Max, Min and sweA.mask
  sweA.mask4[[i]] <- mask(sweA.mask3[[i]], obsMask[[i]])
  sweMax.mask2[[i]] <- mask(sweMax.mask[[i]] , obsMask[[i]])
  sweMin.mask2[[i]] <- mask(sweMin.mask[[i]] , obsMask[[i]])
}



##### Start of melt calculation

#get last day of within 80 % of max
thrsh80 <- function(x, y){
  ifelse(x >= 0.8*y,1,NA )
}

swe80 <- list()
swe80Table <- list()
sweADoy <- list()
for(i in 1:NYears){
  #get function if swe is within the 80% threshold
  swe80[[i]] <- overlay(sweA.mask4[[i]],sweMax.mask2[[i]], fun=thrsh80)
  #put in table
  swe80Table[[i]] <- getValues(swe80[[i]])
  #get doy for swe in data set
  sweADoy[[i]] <- yday(as.Date(names(swePeriod[[i]]), "X%Y.%m.%d"))

}

#function to get last day in column equal to 1
layer80 <- function(x){
  max(which(x == 1))
}


swe80Layer <- list()
swe80LayerA <- list()
sweValues <- list()
startLayers <- list()
sweVs <- numeric()
sweStart  <- list()
for(i in 1:NYears){
  #apply to table
  swe80Layer[[i]] <- apply(swe80Table[[i]],1,layer80)
  #max function returns -INF if NA, fix to NA
  swe80LayerA[[i]] <- ifelse(swe80Layer[[i]] == -Inf,NA,swe80Layer[[i]])
  #get layer with swe
  sweValues[[i]] <- getValues(sweA.mask4[[i]])
  #get swe from layer of last day
  #get first NA if all missing
  startLayers[[i]] <- ifelse(is.na(swe80LayerA[[i]]),1,swe80LayerA[[i]])
  
  sweVs <- numeric()
  for(k in 1:nrow(sweValues[[i]])){
    sweVs[k] <- sweValues[[i]][k,startLayers[[i]][k]]
  }

  #add back into raster
  sweStart[[i]] <- setValues(sweMax.mask2[[i]],sweVs )
  
}


#get layer day
meltStartDay <- numeric()

meltStart  <- list()

for(i in 1:NYears){

    meltStartDay <- ifelse(is.na(swe80LayerA[[i]]),NA,sweADoy[[i]][swe80LayerA[[i]]])


  #add back into the raster
  meltStart[[i]] <- setValues(sweMax.mask2[[i]],meltStartDay )
}  


##### End of melt calculation



#get first day within 20 % of max
thrsh20 <- function(x, y){
  ifelse(x <= 0.2*y,1,NA )
}

swe20 <- list()
swe20Table <- list()
for(i in 1:NYears){
  #get function if swe is within the 20% threshold
  swe20[[i]] <- overlay(sweA.mask4[[i]],sweMax.mask2[[i]], fun=thrsh20)
  #put in table
  swe20Table[[i]] <- getValues(swe20[[i]])
}

#function to get first day in column equal to 1
layer20 <- function(x){
 min(which(x == 1))
}

swe20Layer <- list()
swe20LayerA <- list()
endLayers <- list()
sweEnd <- list()
sweVe <- numeric()
for(i in 1:NYears){
  #apply to table
  swe20Layer[[i]] <- apply(swe20Table[[i]],1,layer20)
  #min will return +Inf for all NA
  swe20LayerA[[i]] <- ifelse(swe20Layer[[i]] == Inf,NA,swe20Layer[[i]])
  
  #get swe from layer of last day
  #get first NA if all missing
  endLayers[[i]] <- ifelse(is.na(swe20LayerA[[i]]),1,swe20LayerA[[i]])
  sweVe <- numeric()
  for(k in 1:nrow(sweValues[[i]])){
    sweVe[k] <- sweValues[[i]][k,endLayers[[i]][k]]
  }
  #add back into raster
  sweEnd[[i]] <- setValues(sweMax.mask2[[i]],sweVe )

}


#get layer day
meltEndDay <- numeric()
meltEnd <- list()

for(i in 1:NYears){
  
  meltEndDay <- ifelse(is.na(swe20LayerA[[i]]),NA,sweADoy[[i]][swe20LayerA[[i]]])
  #add back into the raster
  meltEnd[[i]] <- setValues(sweMax.mask2[[i]],meltEndDay )
}


##### Melt rate calculation


MeltPeriod <- list()
MeltPeriodm <- list()
sweDecline <- list()
sweDeclinem <- list()
Melt.m.day <- list()
Melt.mm.day<- list()

for(i in 1:NYears){  
#duration of melt
MeltPeriod[[i]] <- meltEnd[[i]] - meltStart[[i]]
###### Filter point: exclude melt period with less than 5 days
#unreliable for this analysis and outlier

MeltPeriodm[[i]] <- reclassify(MeltPeriod[[i]], c(-Inf,5,NA))
#total swe decline
sweDecline[[i]] <- sweEnd[[i]] - sweStart[[i]]
#make sure that any melt rates of zero or increases in snow would be excluded
sweDeclinem[[i]] <- reclassify(sweDecline[[i]] , c(0,Inf,NA), include.lowest=TRUE)
}



Melt.m.day <- list()
Melt.mm.day<- list()
for(i in 1:NYears){  
#melt rate in meter per day
Melt.m.day[[i]] <- sweDeclinem[[i]]/MeltPeriodm[[i]]
#melt rate in millimeter per day
Melt.mm.day[[i]] <- (sweDeclinem[[i]]*1000)/MeltPeriodm[[i]]
}

##organize output
#daily swe in EASE projection All
dailySwe <- sweA.ease
#raw swe with mask
dailySwe.mask <- sweA.mask4 

#stacks out melt and calc
melt.mm.day <- stack(Melt.mm.day)
names(melt.mm.day ) <- paste("year",seq(2000,2009))


#length of melt
meltDuration <- stack(MeltPeriodm)
names(meltDuration ) <- paste("year",seq(2000,2009))
#doy start if melt
doyStart<- stack(meltStart)
names(doyStart ) <- paste("year",seq(2000,2009))
#glc info
glc2000 <- glc.maj2
glcID <- data.frame(glc=c(4,5,6,12,13),
                    c("Tree Cover, needle-leaved, evergreen",
                      "Tree Cover, needle-leaved, deciduous",
                      "Tree Cover, mixed leaf type",
                      "Shrub Cover, closed-open, deciduous",
                      "Herbaceous Cover, closed-open"))
#swe Max
maxSwe <- stack( sweMax.mask2)
names(maxSwe) <- paste("year",seq(2000,2009))
plot(maxSwe)

#organize data frame for analysis
#need vcf and air temp
dataAll2000 <- stack(melt.mm.day[[1]],glc2000,doyStart[[1]],maxSwe[[1]])
names(dataAll2000) <- c("melt.mm.day","glc","doyStart","maxSwe.m")
plot(dataAll2000)
dataDF <- getValues(dataAll2000)
coordinatesDF <- data.frame(coordinates(glc2000) )
coordinatesSP <- SpatialPoints(unique(data.frame(x=coordinatesDF$x,y=coordinatesDF$y)), CRS(laea))
#transform to lat/long
coordinatesLL <- spTransform(coordinatesSP, "+init=epsg:4326")
plot(coordinatesLL)

LatLong <- data.frame(lat=coordinatesLL@coords[,2],lon=coordinatesLL@coords[,1])

dataAll2000b <- cbind(dataDF,LatLong)
dataAll2000c <- na.omit(dataAll2000b)
plot(dataAll2000c$lon,dataAll2000c$lat)
#convert to lat long

dataDFc <- cbind(dataDF, coordinatesDF)

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

rm(list=setdiff(ls(), c("dailySwe",
                        "dailySwe.mask",
                        "melt.mm.day",
                        "meltDuration",
                        "doyStart",
                        "glc2000",
                        "glcID",
                        "maxSwe")))


test1 <- stack(Melt.mm.day)
names(test1) <- paste("year",seq(2000,2009))
plot(test1)

plot(sweDecline[[1]])

plot(sweDecline)

tm_shape(Melt.cm.day)+
  tm_raster(title="2000 melt rate (cm/day)", palette= "BuPu", style = "quantile")+
  tm_layout(legend.outside=TRUE)

tm_shape(MeltPeriod)+
  tm_raster(title="2000 melt days", palette= "BuPu",style="quantile")+
  tm_layout(legend.outside=TRUE)

tm_shape(sweMax.mask)+
  tm_raster(title="2000 Maximum swe", palette= "BuPu",style="quantile")+
  tm_layout(legend.outside=TRUE)
#create a mask for snow extent

#raw swe reprojected
sweA.ease
#raw swe with mask
sweA.mask4
#swe at start of melt
sweStart
#swe at the end of melt
sweEnd
#doy of melt end
meltEnd 
#doy of mel start
meltStart
