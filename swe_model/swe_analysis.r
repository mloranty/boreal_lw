##########################################################
#### Analysis of 9 years of swe depletion in boreal   ####
#### forests. Model parameters describing swe         ####
#### depletion are read in. There are three main      ####
#### dataframes: datSwe: the actual swe data          ####
#### b0Out: the rate of swe decline and midOut: the   ####
#### timing of half swe,muB0Out=glc mean slope        ####
#### muMidOut= glc mean midpoint                      ####
##########################################################


###############################################
### read in swe depletion curve output      ###
###############################################


source("c:\\Users\\hkropp\\Documents\\GitHub\\boreal_lw\\swe_model\\swe_output_process.r")
###############################################
### libraries                               ###
###############################################
library(raster)
library(ncdf4)
library(rgdal)
library(gdalUtils)
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"


###############################################
### set up a dataframe with all of the      ###
### data and parameters by cell             ###
###############################################
sweCell <- unique(data.frame(cell=datSwe$cell,gcID=datSwe$gcID,pixID=datSwe$pixID,year=datSwe$year,
								lat=datSwe$y.coord,lon=datSwe$x.coord,
								vcf=datSwe$vcf,zone=datSwe$zone,dayMax=datSwe$dayMax))

#calculate the average air temp during the melt period
tempMelt <- aggregate(datSwe$t.air, by=list(datSwe$gcID,datSwe$pixID,datSwe$year), FUN="mean")
colnames(tempMelt) <- c("gcID","pixID","year","tempK")
tempMelt$tempC <- tempMelt$tempK-273.15
#combine back into sweCell
sweCell2 <- join(sweCell,tempMelt, by=c("gcID","pixID","year"),type="left")

#now join to each output
b0All <- join(b0Out,sweCell2,by=c("year","gcID","pixID"),type="inner")
midAll <- join(midOut,sweCell2,by=c("year","gcID","pixID"),type="inner")

#organize output by year
yearDF <- data.frame(year=unique(sweCell$year))
bOutL <- list()
midOutL <- list()
for(i in 1:dim(yearDF)[1]){
	bOutL[[i]] <- b0All[b0All$year==yearDF$year[i],]
	midOutL[[i]] <- midAll[midAll$year==yearDF$year[i],]
}


###############################################
### set up information for mapping          ###
###############################################
# define the projection - EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 

swe.files <- list.files(pattern =".nc",path =paste0(swepath,"\\swe_mudryk_blended"),full.names=T)

# read one file in to use for reprojecting
pr <- raster(swe.files[1])
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')


## read in file from first year
swe <- raster(swe.files[3])

# crop to boreal region
swe <- crop(swe,c(-180,180,50,90))
# reproject
swe <- projectRaster(swe,pr)

#get the cell to match up to
sweCells <- ncell(swe)
sweCellDF <- data.frame(cell=seq(1,sweCells))

#join back to the swe cell id allowing others to turn to NA
bSwe <- list()
midSwe <- list()
for(i in 1:dim(yearDF)[1]){
	bSwe[[i]] <- join(sweCellDF,bOutL[[i]], by="cell",type="left")
	midSwe[[i]] <- join(sweCellDF,midOutL[[i]], by="cell",type="left")
}
###############################################
### map results                             ###
###############################################
mid2000 <- setValues(swe,midSwe[[1]]$Mean)
mid2001 <- setValues(swe,midSwe[[2]]$Mean)
mid2002 <- setValues(swe,midSwe[[3]]$Mean)
mid2003 <- setValues(swe,midSwe[[4]]$Mean)
mid2004 <- setValues(swe,midSwe[[5]]$Mean)
mid2005 <- setValues(swe,midSwe[[6]]$Mean)
mid2006 <- setValues(swe,midSwe[[7]]$Mean)
mid2007 <- setValues(swe,midSwe[[8]]$Mean)
mid2008 <- setValues(swe,midSwe[[9]]$Mean)
mid2009 <- setValues(swe,midSwe[[10]]$Mean)

b2000 <- setValues(swe,bSwe[[1]]$Mean)
b2001 <- setValues(swe,bSwe[[2]]$Mean)
b2002 <- setValues(swe,bSwe[[3]]$Mean)
b2003 <- setValues(swe,bSwe[[4]]$Mean)
b2004 <- setValues(swe,bSwe[[5]]$Mean)
b2005 <- setValues(swe,bSwe[[6]]$Mean)
b2006 <- setValues(swe,bSwe[[7]]$Mean)
b2007 <- setValues(swe,bSwe[[8]]$Mean)
b2008 <- setValues(swe,bSwe[[9]]$Mean)
b2009 <- setValues(swe,bSwe[[10]]$Mean)

plot(mid2000)
plot(mid2001)
plot(mid2002)
plot(mid2003)
plot(mid2004)
plot(mid2005)
plot(mid2006)
plot(mid2007)
plot(mid2008)
plot(mid2009)

plot(b2000)
plot(b2001)
plot(b2002)
plot(b2003)
plot(b2004)
plot(b2005)
plot(b2006)
plot(b2007)
plot(b2008)
plot(b2009)


#plot each years curve on the map






#plot the means for each glc




#jags regression
plot(midAll$lat,midAll$Mean)
