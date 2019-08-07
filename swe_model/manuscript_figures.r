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
library(maps)
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

modDir <- "z:\\projects\\boreal_swe_depletion\\analysis\\run8"
plotDI <- "c:\\Users\\hkropp\\Google Drive\\Picker17\\figures"

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

#set into raster
rasterCanopy <- setValues(swe,MapCanopy$vcf)
rasterVege <- setValues(swe,MapVege$gcID)
rasterMaxMean <- setValues(swe,MapMax$sweMax)
rasterMaxMissing <- setValues(swe,MapMax$sweMaxExc)
#max missing is predominately zero throughout the entire map. Not worth showing
worldmap <- map("world", ylim=c(40,90), fill=TRUE)
#focus on a smaller extent
worldmap2 <- map("world", ylim=c(50,90))

#world map
world <- project(matrix(c(worldmap$x,worldmap$y), ncol=2,byrow=FALSE),laea)
world2 <- project(matrix(c(worldmap2$x,worldmap2$y), ncol=2,byrow=FALSE),laea)


###############################################
### map results                             ###
###############################################
treePallete <- c(rgb(237,248,233,max=255),
				rgb(199,233,192,max=255),
				rgb(161,217,155,max=255),
				rgb(116,196,118,max=255),
				rgb(49,163,84,max=255),
				rgb(0,109,44,max=255))
vegePallete <- c(rgb(170/255,190/255,140/255),	
				rgb(60/255,60/255,110/255),
				rgb(130/255,160/255,190/255),
				rgb(50/255,80/255,10/255),
				rgb(250/255,120/255,80/255))
				
swePallete <- c(rgb(237,248,251,max=255),
				rgb(191,211,230,max=255),
				rgb(158,188,218,max=255),
				rgb(140,150,198,max=255),
				rgb(136,86,167,max=255),
				rgb(129,15,124,max=255))

hd <- 10
wd1 <- 10
wd2 <- 3
water <- rgb(149/255,218/255,255/255,.3)
land <- rgb(250,230,190, max=255)

png(paste0(plotDI,"\\data_maps.png"), width = 17, height = 5, units = "in", res=300)
	layout(matrix(seq(1,6),ncol=6), width=c(lcm(wd1),lcm(wd2),lcm(wd1),lcm(wd2),lcm(wd1),lcm(wd2)),height=lcm(hd))
	#set up empty plot
	### plot 1 vegetation type ###
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4100000,4100000),ylim=c(-4100000,4100000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterVege,breaks=c(0,1,2,3,4,5),col=vegePallete,add=TRUE )
	### plot 1 legend ###
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	### plot 2 canopy cover ###
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4100000,4100000),ylim=c(-4100000,4100000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterCanopy, breaks=c(0,5,15,25,35,45,75), col=treePallete, add=TRUE)

	### plot 2 legend ###
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1))
	### plot 3 canopy cover ###
		par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4100000,4100000),ylim=c(-4100000,4100000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterMaxMean,breaks=c(0,0.05,0.1,0.2,0.3,0.4,0.5), col=swePallete, add=TRUE)
	### plot 2 legend ###
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1))
dev.off()
