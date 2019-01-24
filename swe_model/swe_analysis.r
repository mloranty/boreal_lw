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
b0All <- join(b0Out,sweCell2,by=c("year","gcID","pixID"),type="left")
midAll <- join(midOut,sweCell2,by=c("year","gcID","pixID"),type="left")

###############################################
### set up information for mapping           ###
###############################################
# define the projection - EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
swe.files <- list.files(pattern =".nc",path ="\\swe_mudryk_blended\\",full.names=T)

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

#need to make a data frame for each year

#join back to the swe cell id allowing others to turn to NA









#plot each years curve on the map






#plot the means for each glc













#jags regression
