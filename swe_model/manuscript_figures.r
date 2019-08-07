##########################################################
#### Analysis of 9 years of swe depletion in boreal   ####
#### forests. Model parameters describing swe         ####
#### depletion are read in. There are three main      ####
#### dataframes: datSwe: the actual swe data          ####
#### b0Out: the rate of swe decline and midOut: the   ####
#### timing of half swe,muB0Out=glc mean slope        ####
#### muMidOut= glc mean midpoint                      ####
#### halfOut= #of days between last Max day and mid   ####
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
library(sp)
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

modDir <- "z:\\projects\\boreal_swe_depletion\\analysis\\run8"


###############################################
### organize all data                       ###
###############################################
sweCell <- unique(data.frame(cell=datSwe$cell,gcID=datSwe$gcID,pixID=datSwe$pixID,year=datSwe$year,
								y=datSwe$y.coord,x=datSwe$x.coord,
								vcf=datSwe$vcf,zone=datSwe$zone,dayMax=datSwe$dayMax, newpixID=datSwe$newpixID))
								
								
#calculate the melt period		
#join midpoint into swe cell
colnames(midOut)[1:7] <-  paste0(colnames(midOut)[1:7],"M")
sweCell2 <- join(sweCell,midOut,by=c("newpixID","gcID","year"),type="left")
colnames(halfOut)[1:7] <- paste0(colnames(halfOut)[1:7],"H")
sweCell3 <- join(sweCell2,halfOut,by=c("newpixID","gcID","year"),type="left")

#### data filter #####
#there are some very fast melt periods
#where the standard deviation of the midpoint
#is within the range of the onset. Looking at the
#melt period won't be reliable. 					
sweCell3 <- sweCell3[round(sweCell3$MeanM-sweCell3$SDM)>sweCell3$dayMax,]

#parm all
parmAll <- join(b0Out,sweCell3,by=c("year","gcID","newpixID"),type="inner")


################################################################################
################################################################################
############### Figure 1. Map of data inputs                     ############### 
###############including % tree cover, vegetation class          ###############
############### and average maximum swe during period            ###############
################################################################################
################################################################################
###############################################
### organize descriptive data               ###
###############################################

#aggregate swe max, but first need to pull it out from join to full data
sweMaxDF <- unique(data.frame(cell=sweAll$cell, year=sweAll$year,sweMax=sweAll$sweMax))
#aggregate by each cell
sweMaxSumm <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="mean")
colnames(sweMaxSumm) <- c("cell","sweMax")
sweMaxSumm$sweMaxSD <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="sd")$x
sweMaxSumm$sweMaxN <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="length")$x
sweMaxSumm$sweMaxExc <- 10 - sweMaxSumm$sweMaxN 
#get canopy cover and gcID 

canopyCov <- unique(data.frame(cell=sweAll$cell, vcf=sweAll$vcf))
#organize vege type
gcIDSumm <- unique(data.frame(cell=sweAll$cell, gcID=sweAll$gcID))
#get names for plot
gcNames <- unique(data.frame(names=sweAll$names, gcID=sweAll$gcID))
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
MapCanopy <- join(sweCellDF,canopyCov, by="cell",type="left")
MapVege <- join(sweCellDF,gcIDSumm, by="cell",type="left")
MapMax <- join(sweCellDF,sweMaxSumm, by="cell",type="left")
###############################################
### map results                             ###
###############################################


rasterCanopy <- setValues(swe,MapCanopy$vcf)
rasterVege <- setValues(swe,MapVege$gcID)
rasterMaxMean <- setValues(swe,MapMax$sweMax)
rasterMaxMissing <- setValues(swe,MapMax$sweMaxExc)


