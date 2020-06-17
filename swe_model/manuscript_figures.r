##########################################################
#### Analysis of 9 years of swe depletion in boreal   ####
#### forests.                                         ####
##########################################################


###############################################
### read in swe depletion curve output      ###
###############################################


source("c:\\Users\\hkropp\\Documents\\GitHub\\boreal_lw\\swe_model\\swe_data_org.r")

###############################################
### libraries                               ###
###############################################
library(raster)
library(ncdf4)
library(rgdal)
library(gdalUtils)
library(sp)
library(maps)
library(rgeos)
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

modDir <- "z:\\projects\\boreal_swe_depletion\\analysis\\run19"
plotDI <- "c:\\Users\\hkropp\\Google Drive\\Picker17\\figures"

###############################################
### define color palette                    ###
###############################################
vegePallete <- c(rgb(50/255,80/255,10/255), #evergreen needleleaf,	
				rgb(250/255,120/255,80/255), #mixed boreal,
				rgb(170/255,190/255,140/255), #herbaceous,
				rgb(130/255,160/255,190/255),# deciduous needleleaf,
				rgb(60/255,60/255,110/255)) #deciduous shrub)
vegePallete2 <-	c(rgb(50/255,80/255,10/255,.1),	
				rgb(250/255,120/255,80/255,.1),
				rgb(170/255,190/255,140/255,.1),
				rgb(130/255,160/255,190/255,.1),
				rgb(60/255,60/255,110/255,.1))		
vegePallete3 <-	c(rgb(50/255,80/255,10/255,.5),	
				rgb(250/255,120/255,80/255,.5),
				rgb(170/255,190/255,140/255,.5),
				rgb(130/255,160/255,190/255,.5),
				rgb(60/255,60/255,110/255,.5))		

#add better names
IDSglc$name2 <- c("Evergreen needleleaf", "Mixed boreal", "Herbaceous","Deciduous needleleaf","Deciduous shrub")

				
###############################################
### read in regression output               ###
###############################################
#read in model output
datS <- read.csv(paste0(modDir,"\\curve_mod_stats.csv"))
datQ <- read.csv(paste0(modDir,"\\curve_mod_quant.csv"))

#combine data frames
datC <- cbind(datS,datQ)
#pull out parameter names
dexps<-"\\[*[[:digit:]]*\\]"
datC$parm <- gsub(dexps,"",rownames(datC))
unique(datC$parm)
datC$parm2 <- gsub("\\d","",rownames(datC))

#pull out parameters
#transformed intercepts
beta0NL <- datC[datC$parm == "trB0",] 
#nontransformed regression parameters
beta0 <- datC[datC$parm == "betaB0S",] 
beta1 <- datC[datC$parm == "betaB1",] 
beta2 <- datC[datC$parm == "betaB2",] 
beta3 <- datC[datC$parm == "betaB3",] 
beta4 <- datC[datC$parm == "betaB4",] 
#add indicator if parameter is significant
beta1$sig <- ifelse(beta1$X2.5.<0&beta1$X97.5.<0,1,
						ifelse(beta1$X2.5.>0&beta1$X97.5.>0,1,0))
						
						
beta3$sig <- ifelse(beta3$X2.5.<0&beta3$X97.5.<0,1,
						ifelse(beta3$X2.5.>0&beta3$X97.5.>0,1,0))					
beta4$sig <- ifelse(beta4$X2.5.<0&beta4$X97.5.<0,1,
						ifelse(beta4$X2.5.>0&beta4$X97.5.>0,1,0))	
						
#check if any nonsignificant slopes to account for

length(which(beta1$sig == 0))						
						
length(which(beta3$sig == 0))	
	

#pull out regression means for plotting

mu.Temp <- datC[datC$parm2 == "mu.Temp[,]",]
mu.Onset <- datC[datC$parm2 == "mu.Onset[,]",]
mu.Lat <- datC[datC$parm2 == "mu.Lat[,]",]
						
###############################################
### add in unique id for model              ###
###############################################

#join unique glc id gcID
cellSwe2 <- join(cellSwe, IDSglc, by="zone",type="left")


#create unique year ID
IDSyears <- unique(data.frame(year=cellSwe2$year))								
IDSyears$yearID <- seq(1,nrow(IDSyears))
#join back into cellSwe
cellSwe3 <- join(cellSwe2, IDSyears, by="year",type="left")

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
meltSwe <- list()

for(i in 1:dim(IDSyears)[1]){
	meltSwe[[i]] <- join(sweCellDF,cellSwe3[cellSwe3$year==IDSyears$year[i],], by="cell",type="left")

}
###############################################
### map results                             ###
###############################################


#get lat long for each cell
#create a spatial points
sweSP <- SpatialPoints(unique(data.frame(x=sweAll$x.coord,y=sweAll$y.coord,cell=sweAll$cell)), CRS(laea))
#transform for wgs lat long
sweSPr <- spTransform(sweSP, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

sweLL <- data.frame(sweSPr@coords)
colnames(sweLL) <- c("Lon","Lat","cell")
#join with xy coord
sweSpatial <- unique(data.frame(x=sweAll$x.coord,y=sweAll$y.coord,cell=sweAll$cell))

#join back into dataframe of results
cellSwe4 <- join(cellSwe3,sweLL, by="cell",type="left")
cellSwe5 <- join(cellSwe4,sweSpatial, by="cell",type="left")
###############################################
### finish organizing model                 ###
###############################################
#need to organize table for eps ids
epsTable <- unique(data.frame(gcID=cellSwe5$gcID,year=cellSwe5$year))
epsTable <- epsTable[order(epsTable$gcID,epsTable$year),]
#this order will be by GCID
epsTable$gcyearID <- seq(1,dim(epsTable)[1])

#join back into b0
cellSwe6 <- join(cellSwe5,epsTable, by=c("gcID","year"),type="left")

#create index for averaging eps
gcIndT <- unique(data.frame(gcID=epsTable$gcID))
startID <- numeric(0)
endID <- numeric(0)

for(i in 1:dim(gcIndT)[1]){
		startID[i] <- head(which(epsTable$gcID==gcIndT$gcID[i]))[1]
		endID [i] <- tail(which(epsTable$gcID==gcIndT$gcID[i]))[6]
}

#index for spatial random effect
cellDF <- unique(data.frame(cell=cellSwe6$cell, y=cellSwe6$y,x=cellSwe6$x,gcID=cellSwe6$gcID))
cellDF$cellID <- seq(1,nrow(cellDF))

#join back into b0All
cellSwe7 <- join(cellSwe6, cellDF, by=c("cell","x","y","gcID"), type="left")

cellSwe7$absRate <- abs(cellSwe7$meltRateCM)
cellSwe7$logAbsRate <- log(cellSwe7$absRate )
sweRate <- cellSwe7


#set up data for plotting
tempMean <- seq(floor(range(cellSwe7$tair)[1]),ceiling(range(cellSwe7$tair)[2]), length.out=200)
CanopyMean <- seq(floor(range(cellSwe7$vcf)[1]),ceiling(range(cellSwe7$vcf)[2]), length.out=200)
SdayMean <- seq(floor(range(cellSwe7$meltStart)[1]),ceiling(range(cellSwe7$meltStart)[2]), length.out=200)
LatMean <- seq(50,75, length.out=200)
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

#get average swe max for each of the cells
sweMaxDF <- unique(data.frame(cell=sweRate$cell, year=sweRate$year,sweMax=sweRate$sweMax))
#aggregate by each cell
sweMaxSumm <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="mean")
colnames(sweMaxSumm) <- c("cell","sweMax")
sweMaxSumm$sweMaxSD <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="sd")$x
sweMaxSumm$sweMaxN <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="length")$x

canopyCov <- unique(data.frame(cell=sweRate$cell, vcf=sweRate$vcf))
#organize vege type
gcIDSumm <- unique(data.frame(cell=sweRate$cell, gcID=sweRate$gcID))
#get names for plot
gcNames <- unique(data.frame(names=sweRate$names, gcID=sweRate$gcID))


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

#max missing is predominately zero throughout the entire map. Not worth showing
worldmap <- map("world", ylim=c(40,90), fill=TRUE)
#focus on a smaller extent
worldmap2 <- map("world", ylim=c(50,90))

#world map
world <- project(matrix(c(worldmap$x,worldmap$y), ncol=2,byrow=FALSE),laea)
world2 <- project(matrix(c(worldmap2$x,worldmap2$y), ncol=2,byrow=FALSE),laea)

###make polygon to cover up non study area####
#make a point at center of map
pt1 <- SpatialPoints(data.frame(x=0,y=0), CRS(laea))
#make a buffer around center of plot choosing distance that is relevant
ptBuff <- buffer(pt1,4500000)
#set up bounds to extend beyond study area
xcor <- c(-8000000,-8000000,8000000,8000000)
ycor <- c(-8000000,8000000,8000000,-8000000)
#make empty plot polygon
boxC <- cbind(xcor,ycor)
p <- Polygon(boxC)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
proj4string(sps) = CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" )
#remove study area from empty plot
PolyBlock <- gDifference(sps,ptBuff, byid=TRUE)


###############################################
### map results                             ###
###############################################


treePallete <- c(rgb(229,245,224,max=255),
				rgb(199,233,192,max=255),
				rgb(161,217,155,max=255),
				rgb(116,196,118,max=255),
				rgb(65,171,93,max=255),
				rgb(35,139,69,max=255),
				rgb(0,109,44,max=255),
				rgb(0,68,27,max=255))

				

hd <- 12
wd1 <- 12
wd2 <- 8
water <- rgb(149/255,218/255,255/255,.3)
land <- rgb(250,230,190, max=255)
#size of panel label
mx <- 2
#line for panel label
pll <- .5
#vege type breaks
vegeBr <- c(0,1,2,3,4,5)
canopyBr <- c(0,10,20,30,40,50,60,70,80)
sweBr <-c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5)
#size of axis
cxa <- 1.75

png(paste0(plotDI,"\\data_maps.png"), width = 17, height = 7, units = "in", res=300)
	layout(matrix(seq(1,4),ncol=4), width=c(lcm(wd1),lcm(wd2),lcm(wd1),lcm(wd2)),height=lcm(hd))
	#set up empty plot
	### plot 1 vegetation type ###
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterVege,breaks=vegeBr,col=vegePallete,add=TRUE )
	mtext("A",at=4100000,side=2,line=pll, las=2,cex=mx)
	plot(PolyBlock, col="white",border="white", add=TRUE)
	### plot 1 legend ###
	par(mai=c(0,.25,0,2))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(vegeBr)-1)){
		polygon(c(0,0,1,1), 
		c(vegeBr[i]/(length(vegeBr)-1),vegeBr[i+1]/(length(vegeBr)-1),vegeBr[i+1]/(length(vegeBr)-1),vegeBr[i]/(length(vegeBr)-1)),
		col=vegePallete[i],border=NA)
	}
	axis(4,(vegeBr[1:5]/5)+.1,IDSglc$name2,cex.axis=cxa,las=2)

	
	### plot 2 canopy cover ###
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterCanopy, breaks=canopyBr, col=treePallete, add=TRUE)
	mtext("B",at=4100000,side=2,line=pll, las=2,cex=mx)
	plot(PolyBlock, col="white",border="white", add=TRUE)
	### plot 2 legend ###
	par(mai=c(0,.25,0,2))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(canopyBr)-1)){
		polygon(c(0,0,1,1), 
			c(canopyBr[i]/canopyBr[length(canopyBr)],canopyBr[i+1]/canopyBr[length(canopyBr)],canopyBr[i+1]/canopyBr[length(canopyBr)],canopyBr[i]/canopyBr[length(canopyBr)]),
			col=treePallete[i],border=NA)
	}
	axis(4,canopyBr/canopyBr[length(canopyBr)],canopyBr,cex.axis=cxa,las=2)	
	
dev.off()





################################################################################
################################################################################
############### Figure 2. Map of data related to melt            ############### 
################################################################################
################################################################################
#get average, sd of swe rate by cell, then show average maximum and average onset


###############################################
### organize descriptive data               ###
###############################################

#get average swe max for each of the cells
sweMaxDF <- unique(data.frame(cell=sweRate$cell, year=sweRate$year,sweMax=sweRate$sweMax))
#aggregate by each cell
sweMaxSumm <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="mean")
colnames(sweMaxSumm) <- c("cell","sweMax")
sweMaxSumm$sweMaxSD <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="sd")$x
sweMaxSumm$sweMaxN <- aggregate(sweMaxDF$sweMax, by=list(sweMaxDF$cell), FUN="length")$x

#get average swe max for each of the cells
onsetDF <- unique(data.frame(cell=sweRate$cell, year=sweRate$year,sweOnset=sweRate$meltStart))
#aggregate by each cell
onsetSumm <- aggregate(onsetDF$sweOnset, by=list(onsetDF$cell), FUN="mean")
colnames(onsetSumm) <- c("cell","sweOnset")

#get average melt rate
sweMeltSumm <- aggregate(sweRate$absRate, by=list(sweRate$cell), FUN="mean")
colnames(sweMeltSumm) <- c("cell","aveMelt")
sweMeltSumm$CV <- aggregate(sweRate$absRate, by=list(sweRate$cell), FUN="sd")$x/sweMeltSumm$aveMelt

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

MapMax <- join(sweCellDF,sweMaxSumm, by="cell",type="left")
MapMelt <- join(sweCellDF,sweMeltSumm, by="cell",type="left")
MapOnset <- join(sweCellDF,onsetSumm, by="cell",type="left")


#set into raster

rasterMaxMean <- setValues(swe,MapMax$sweMax)
rasterMeltMean <- setValues(swe,MapMelt$aveMelt)
rasterMeltCV <- setValues(swe,MapMelt$CV)
rasterOnset <- setValues(swe,MapOnset$sweOnset)

#max missing is predominately zero throughout the entire map. Not worth showing
worldmap <- map("world", ylim=c(40,90), fill=TRUE)
#focus on a smaller extent
worldmap2 <- map("world", ylim=c(50,90))

#world map
world <- project(matrix(c(worldmap$x,worldmap$y), ncol=2,byrow=FALSE),laea)
world2 <- project(matrix(c(worldmap2$x,worldmap2$y), ncol=2,byrow=FALSE),laea)

###make polygon to cover up non study area####
#make a point at center of map
pt1 <- SpatialPoints(data.frame(x=0,y=0), CRS(laea))
#make a buffer around center of plot choosing distance that is relevant
ptBuff <- buffer(pt1,4500000)
#set up bounds to extend beyond study area
xcor <- c(-8000000,-8000000,8000000,8000000)
ycor <- c(-8000000,8000000,8000000,-8000000)
#make empty plot polygon
boxC <- cbind(xcor,ycor)
p <- Polygon(boxC)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
proj4string(sps) = CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" )
#remove study area from empty plot
PolyBlock <- gDifference(sps,ptBuff, byid=TRUE)
###############################################
### map results                             ###
###############################################


swemaxPallete <- c(rgb(247,252,253,max=255),
				rgb(224,236,244,max=255),
				rgb(158,188,218,max=255),
				rgb(140,150,198,max=255),
				rgb(140,107,177,max=255),
				rgb(136,65,157,max=255),
				rgb(129,15,124,max=255),
				rgb(77,0,75,max=255))
				
SDPallete <- c(	rgb(255,245,240,max=255),
				rgb(254,224,210,max=255),
				rgb(252,187,161,max=255),
				rgb(252,146,114,max=255),
				rgb(251,106,74,max=255),
				rgb(239,59,44,max=255),
				rgb(203,24,29,max=255),
				rgb(153,0,13,max=255))
				
OnsetPallete <- c(rgb(255,247,251,max=255),
					rgb(236,226,240,max=255),
					rgb(208,209,230,max=255),
					rgb(166,189,219,max=255),
					rgb(103,169,207,max=255),
					rgb(54,144,192,max=255),
					rgb(2,129,138,max=255),
					rgb(1,100,80,max=255))

swePallete <- rev(c(rgb(178,24,43,max=255),
				rgb(214,96,77,max=255),
				rgb(244,165,130,max=255),
				rgb(253,219,199,max=255),
				rgb(209,229,240,max=255),
				rgb(146,197,222,max=255),
				rgb(67,147,195,max=255),
				rgb(33,102,172,max=255)))
				
sweSDBr <- round(seq(0,0.95, length.out=9),2)
sweMaxBr <- c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5)				
OnsetBr <- round(seq(45,155, length.out=9))
sweBr <- round(seq(0.04,2, length.out=9),2) 

				
hd <- 14
wd1 <- 14
wd2 <- 8
water <- rgb(149/255,218/255,255/255,.3)
land <- rgb(250,230,190, max=255)
#size of panel label
mx <- 2
#line for panel label
pll <- .5
#size of axis
cxa <- 2			
				
png(paste0(plotDI,"\\maps_swe.png"), width = 18, height = 12, units = "in", res=300)
	layout(matrix(seq(1,8),ncol=4, byrow=TRUE), width=c(lcm(wd1),lcm(wd2),lcm(wd1),lcm(wd2)),height=c(lcm(hd),lcm(hd)))	

	### plot 1 swe ave ###
	par(mai=c(.25,.25,.25,.25))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterMeltMean,breaks=sweBr, col=swePallete, add=TRUE)
	mtext("A",at=4100000,side=2,line=pll, las=2,cex=mx)
		plot(PolyBlock, col="white",border="white", add=TRUE)
	### legent plot 1 swe ave ###
	par(mai=c(0.25,0.25,0.25,2))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(sweBr)-1)){
		polygon(c(0,0,1,1), 
			c(sweBr[i]/sweBr[length(sweBr)],
			sweBr[i+1]/sweBr[length(sweBr)],
			sweBr[i+1]/sweBr[length(sweBr)],
			sweBr[i]/sweBr[length(sweBr)]),
			col=swePallete[i],border=NA)
	}
	axis(4,sweBr/sweBr[length(sweBr)],sweBr,cex.axis=cxa,las=2)	
	
	
	### plot 2 swe sd ###
	par(mai=c(.25,.25,.25,.25))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterMeltCV,breaks=sweSDBr, col=SDPallete, add=TRUE)
	mtext("B",at=4100000,side=2,line=pll, las=2,cex=mx)
		plot(PolyBlock, col="white",border="white", add=TRUE)
	### legent plot 1 swe ave ###
	par(mai=c(0.25,0.25,0.25,2))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(sweSDBr)-1)){
		polygon(c(0,0,1,1), 
			c(sweSDBr[i]/sweSDBr[length(sweSDBr)],
			sweSDBr[i+1]/sweSDBr[length(sweSDBr)],
			sweSDBr[i+1]/sweSDBr[length(sweSDBr)],sweSDBr[i]/sweSDBr[length(sweSDBr)]),
			col=SDPallete[i],border=NA)
	}
	axis(4,sweSDBr/sweSDBr[length(sweSDBr)],sweSDBr,cex.axis=cxa,las=2)		
	
	
	### plot 3 swe max ###
		par(mai=c(.25,.25,.25,.25))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterMaxMean,breaks=sweMaxBr, col=swemaxPallete, add=TRUE)
	mtext("C",at=4100000,side=2,line=pll, las=2,cex=mx)
		plot(PolyBlock, col="white",border="white", add=TRUE)
	### plot 3 swe max ###
	par(mai=c(0.25,0.25,0.25,2))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(sweMaxBr)-1)){
		polygon(c(0,0,1,1), 
			c(sweMaxBr[i]/sweMaxBr[length(sweMaxBr)],sweMaxBr[i+1]/sweMaxBr[length(sweMaxBr)],sweMaxBr[i+1]/sweMaxBr[length(sweMaxBr)],sweMaxBr[i]/sweMaxBr[length(sweMaxBr)]),
			col=swemaxPallete[i],border=NA)
	}
	axis(4,sweMaxBr/sweMaxBr[length(sweMaxBr)],sweMaxBr,cex.axis=cxa,las=2)	
	
	
### plot 4 onset ###
		par(mai=c(.25,.25,.25,.25))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
	#color background
	polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
	#boundaries
	points(world, type="l", lwd=2, col="grey65")
	#continent color
	polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
	#plot points
	image(rasterOnset,breaks=OnsetBr, col=OnsetPallete, add=TRUE)
	mtext("D",at=4100000,side=2,line=pll, las=2,cex=mx)
		plot(PolyBlock, col="white",border="white", add=TRUE)
	### plot 4 onset ###
	par(mai=c(0.25,0.25,0.25,2))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(OnsetBr)-1)){
		polygon(c(0,0,1,1), 
			c((OnsetBr[i]-OnsetBr[1])/(OnsetBr[length(OnsetBr)]-OnsetBr[1]),
			(OnsetBr[i+1]-OnsetBr[1])/(OnsetBr[length(OnsetBr)]-OnsetBr[1]),
			(OnsetBr[i+1]-OnsetBr[1])/(OnsetBr[length(OnsetBr)]-OnsetBr[1]),
			(OnsetBr[i]-OnsetBr[1])/(OnsetBr[length(OnsetBr)]-OnsetBr[1])),
			col=OnsetPallete[i],border=NA)
	}
	axis(4,(OnsetBr-OnsetBr[1])/(OnsetBr[length(OnsetBr)]-OnsetBr[1]),OnsetBr,cex.axis=cxa,las=2)	
	
dev.off()



################################################################################
################################################################################
############### Figure 2. Plot of regression with melt rate      ############### 
################################################################################
################################################################################
sweRate
#organize model output
mu.Temp$gcID <- rep(seq(1,5),each=200) 
mu.Onset$gcID <- rep(seq(1,5),each=200) 
mu.Lat$gcID <- rep(seq(1,5),each=200) 

#intercepts

plotTree <- c(1,4,2)	
plotTun <- c(3,5)			
				
wd1 <- 11
hd1 <- 11	
yl <- -3.6
yh <- 1.2	
xl1 <- 	floor(range(sweRate$tair)[1])
xh1 <-	ceiling(range(sweRate$tair)[2])
xl2 <- floor(range(sweRate$vcf)[1])
xh2 <-	ceiling(range(sweRate$vcf)[2])
xl3 <- range(sweRate$meltStart)[1] - 1
xh3 <-	range(sweRate$meltStart)[2] + 1
xl4 <- floor(range(sweRate$Lat)[1])
xh4 <- ceiling(range(sweRate$Lat)[2])
#axis labels
xs1 <- seq(xl1,xh1-3, by=3)
xs2 <- seq(xl2,xh2, by= 15)
xs3 <- seq(60,xh3, by= 30)
ys <- seq(yl,yh, by = 1)
xs4 <- seq(xl4,xh4,by=5)

#width of regression line
mlw <- 4
#width of ticks
tlw <- 4
#axis tick label line
tll <- 2
#axis label size
alc <- 2
#plot label text size
plc <- 3
#x label plot line
xpl <- 6
#size of panel letter
ttx <- 4
#legend size
legcex <- 2.5

dlTemp <- numeric(0)
dhTemp <- numeric(0)
dlvcf<- numeric(0)
dhvcf <- numeric(0)
dlOnset <- numeric(0)
dhOnset <- numeric(0)
dlLat <- numeric(0)
dhLat <- numeric(0)
for(i in 1:5){
	dlTemp[i] <- floor(min(sweRate$tair[sweRate$gcID == i]))
	dhTemp[i] <- ceiling(max(sweRate$tair[sweRate$gcID == i]))
	dlvcf[i] <- floor(min(sweRate$vcf[sweRate$gcID == i]))
	dhvcf[i] <- ceiling(max(sweRate$vcf[sweRate$gcID == i]))
	dlOnset[i] <- floor(min(sweRate$meltStart[sweRate$gcID == i]))
	dhOnset[i] <- ceiling(max(sweRate$meltStart[sweRate$gcID == i]))
	dlLat[i] <- floor(min(sweRate$Lat[sweRate$gcID == i]))
	dhLat[i] <- ceiling(max(sweRate$Lat[sweRate$gcID == i]))
	
	
}




png(paste0(plotDI,"\\regression.png"), width = 55, height = 32, units = "cm", res=300)
	layout(matrix(seq(1,8),ncol=4, byrow=TRUE), width=rep(lcm(wd1),4),height=rep(lcm(hd1),2))
	par(mai=c(0,0,0,0))
	#temperature trees
	plot(c(0,1),c(0,1), type="n", xlim=c(xl1,xh1), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
	for(i in plotTree){
		points(	sweRate$tair[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
	}
	for(i in plotTree){	
		polygon(c(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
					rev(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
				c(mu.Temp$X2.5.[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
					rev(mu.Temp$X97.5[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
				border=NA, col=vegePallete3[i])

		points(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
				mu.Temp$Mean[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
				type="l", lwd=mlw, col=vegePallete[i])
	}
	axis(2, ys, rep(" ",length(ys)), lwd.ticks=tlw)
	mtext(ys,at=ys, line=tll, cex=alc, side=2,las=2)
	box(which="plot")
	mtext(expression(paste("log(Melt Rate (cm day"^"-1","))")), side=2, outer=TRUE,line= -5, cex=plc)
	text(xl1+(.05*(xh1-xl1)), yh-(.05*(yh-yl)), "a", cex=ttx)
	legend("bottomright", paste(IDSglc$name2[plotTree]), col=vegePallete[plotTree],cex=legcex, lwd=mlw,lty=1, bty="n")
	#tree cover trees
	par(mai=c(0,0,0,0))	
	plot(c(0,1),c(0,1), type="n", xlim=c(xl2,xh2), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
	for(i in plotTree){
		points(	sweRate$vcf[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
	}
	for(i in plotTree){	
		polygon(c(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
				rev(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
				c(rep(beta0$X2.5.[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
				rep(beta0$X97.5[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]]))),
				border=NA, col=vegePallete3[i])
		
		points(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
				rep(beta0$Mean[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
				type="l", lwd=mlw, lty=3, col=vegePallete[i])		
				
	}
		box(which="plot")
		text(xl2+(.05*(xh2-xl2)), yh-(.05*(yh-yl)), "c", cex=ttx)
	#day of onset trees
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1), type="n", xlim=c(xl3,xh3), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
	for(i in plotTree){
		points(	sweRate$meltStart[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
	}
	for(i in plotTree){	
		polygon(c(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
					rev(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
				c(mu.Onset$X2.5.[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
					rev(mu.Onset$X97.5[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
				border=NA, col=vegePallete3[i])

		points(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
				mu.Onset$Mean[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
				type="l", lwd=mlw, col=vegePallete[i])
	}
		box(which="plot")
		text(xl3+(.05*(xh3-xl3)), yh-(.05*(yh-yl)), "e", cex=ttx)
		
		#latitude trees
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1), type="n", xlim=c(xl4,xh4), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
	for(i in plotTree){
		points(	sweRate$Lat[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
	}
		for(i in plotTree){	
		
			polygon(c(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]],
				rev(LatMean[LatMean >= dlLat[i] & LatMean<= dhLat[i]])),
				c(rep(beta0$X2.5.[i], length(LatMean[LatMean >= dlLat[i] & LatMean<= dhLat[i]])),
				rep(beta0$X97.5[i], length(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]]))),
				border=NA, col=vegePallete3[i])
		
			points(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]],
				rep(beta0$Mean[i], length(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]])),
				type="l", lwd=mlw, lty=3, col=vegePallete[i])		
		
	}
		box(which="plot")
		text(xl4+(.05*(xh4-xl4)), yh-(.05*(yh-yl)), "g", cex=ttx)	
		
		
	#temperature tundra
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1), type="n", xlim=c(xl1,xh1), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
		for(i in plotTun){
		points(	sweRate$tair[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
		}
		for(i in plotTun){
		polygon(c(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
					rev(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
				c(mu.Temp$X2.5.[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
					rev(mu.Temp$X97.5[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
				border=NA, col=vegePallete3[i])

		points(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
				mu.Temp$Mean[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
				type="l", lwd=mlw, col=vegePallete[i])		
	}	
	axis(2, ys, rep(" ",length(ys)), lwd.ticks=tlw)
	mtext(ys,at=ys, line=tll, cex=alc, side=2,las=2)
	axis(1, xs1, rep(" ",length(xs1)), lwd.ticks=tlw)
	mtext(xs1,at=xs1, line=tll, cex=alc, side=1)
	box(which="plot")
	mtext("Temperature (c)", side=1,line= xpl, cex=plc)
	text(xl1+(.05*(xh1-xl1)), yh-(.05*(yh-yl)), "b", cex=ttx)
	legend("bottomright", paste(IDSglc$name2[plotTun]), col=vegePallete[plotTun],cex=legcex, lwd=mlw,lty=1, bty="n")
	#tree cover tundra
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1), type="n", xlim=c(xl2,xh2), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
	for(i in plotTun){
		points(	sweRate$vcf[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
	}
	for(i in plotTun){	
		polygon(c(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
				rev(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
				c(rep(beta0$X2.5.[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
				rep(beta0$X97.5[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]]))),
				border=NA, col=vegePallete3[i])
		
		points(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
				rep(beta0$Mean[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
				type="l", lwd=mlw, lty=3, col=vegePallete[i])		
				
	}
	axis(1, xs2, rep(" ",length(xs2)), lwd.ticks=tlw)
	mtext(xs2,at=xs2, line=tll, cex=alc, side=1)
	box(which="plot")
	mtext("Canopy cover (%)", side=1,line= xpl, cex=plc)
	text(xl2+(.05*(xh2-xl2)), yh-(.05*(yh-yl)), "d", cex=ttx)
	#onset day tundra
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1), type="n", xlim=c(xl3,xh3), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
	for(i in plotTun){
		points(	sweRate$meltStart[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
	}
	for(i in plotTun){	
		polygon(c(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
					rev(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
				c(mu.Onset$X2.5.[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
					rev(mu.Onset$X97.5[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
				border=NA, col=vegePallete3[i])

		points(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
				mu.Onset$Mean[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
				type="l", lwd=mlw, col=vegePallete[i])
	}
	axis(1, xs3, rep(" ",length(xs3)), lwd.ticks=tlw)
	mtext(xs3,at=xs3, line=tll, cex=alc, side=1)
	box(which="plot")
	mtext("Onset day of year", side=1,line= xpl, cex=plc)
	text(xl3+(.05*(xh3-xl3)), yh-(.05*(yh-yl)), "f", cex=ttx)
	#latitude tundra
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1), type="n", xlim=c(xl4,xh4), ylim=c(yl,yh), xaxs="i",yaxs="i",
		xlab= " ", ylab=" ", axes=FALSE)
	for(i in plotTun){
		points(	sweRate$Lat[sweRate$gcID == i],
				sweRate$logAbsRate[sweRate$gcID == i], col=vegePallete2[i], pch=19)
	}
		for(i in plotTun){	

			polygon(c(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]],
				rev(LatMean[LatMean >= dlLat[i] & LatMean<= dhLat[i]])),
				c(rep(beta0$X2.5.[i], length(LatMean[LatMean >= dlLat[i] & LatMean<= dhLat[i]])),
				rep(beta0$X97.5[i], length(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]]))),
				border=NA, col=vegePallete3[i])
		
			points(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]],
				rep(beta0$Mean[i], length(LatMean[LatMean >= dlLat[i] & LatMean <= dhLat[i]])),
				type="l", lwd=mlw, lty=3, col=vegePallete[i])		
				
	}
		box(which="plot")
		text(xl4+(.05*(xh4-xl4)), yh-(.05*(yh-yl)), "h", cex=ttx)	
		
	axis(1, xs4, rep(" ",length(xs4)), lwd.ticks=tlw)
	mtext(xs4,at=xs4, line=tll, cex=alc, side=1)
	mtext("Latitude", side=1,line= xpl, cex=plc)
dev.off()


################################################################################
################################################################################
############### Figure 3. Plot of regression intercepts          ############### 
################################################################################
################################################################################

beta0NL$gcID <- seq(1,5)
#join matching
intercept <- join(beta0NL,IDSglc, by="gcID", type="left")
plotOrder <- c(1,4,2,3,5)

#break up names 
nameSplit1 <- character(0)
nameSplit2 <- character(0)
for(i in 1:5){
	nameSplit1[i] <- strsplit(intercept$name2[i], " ")[[1]][1]
	nameSplit2[i] <- strsplit(intercept$name2[i], " ")[[1]][2]
}
nameSplit2 <- ifelse(is.na(nameSplit2), " ", nameSplit2)

xseq <- seq(1,10, by=2)

wd1 <- 18
hd1 <- 18
#error bar width
eew <- 1
#mean bar width
mlw <- 2
#tick arrow width
tlw <- 2

png(paste0(plotDI,"\\intercepts1.png"), width = 20, height = 20, units = "cm", res=300)
	layout(matrix(c(1),ncol=1, byrow=TRUE), width=lcm(wd1),height=lcm(hd1))
	plot(c(0,1),c(0,1), xlim=c(0,11), ylim=c(.2,0.6), axes=FALSE, type="n", xlab = " ", ylab= " ",
		xaxs="i", yaxs="i")
		
	for(j in 1:5){
		i <- plotOrder[j]
		arrows(xseq[j],intercept$X2.5.[i],xseq[j],intercept$X97.5.[i], code=0, lwd=eew)
		polygon(c(xseq[j]-.5,xseq[j]-.5,xseq[j]+.5,xseq[j]+.5),
				c(intercept$X25.[i],intercept$X75.[i],intercept$X75.[i],intercept$X25.[i]),
				border=NA,col=vegePallete3[i])
		arrows(	xseq[j]-.5,intercept$Mean[i],xseq[j]+.5,	intercept$Mean[i],code=0,lwd=mlw,
				col=vegePallete[i])

	}
	axis(1, xseq, rep(" ",length(xseq)),lwd.ticks=tlw)
	axis(2, seq(0,.6, by=.1),rep(" ",length(seq(0,.6, by=.1))),lwd.ticks=tlw)
	mtext(paste(nameSplit1[plotOrder]),at=xseq,side=1,line=1,cex=1)
	mtext(paste(nameSplit2[plotOrder]),at=xseq,side=1,line=2,cex=1)
	mtext(seq(.2,.6, by=.1), at=seq(.2,.6, by=.1), side=2, las=2, line=1, cex=1)
	mtext(expression(paste("Melt rate (cm day"^"-1",")")), side=2, line=3, cex=1.5)
	mtext("Landcover type", side=1, line=3, cex=1.5)
dev.off()		


################################################################################
################################################################################
############### Figure 3. violin plot of melt rate               ############### 
################################################################################
################################################################################

#get distributions
histL <- list()
densityH <- numeric(0)
maxS <- numeric(0)

for(i in 1:5){

	histL[[i]] <- hist(sweRate$absRate[sweRate$gcID == i], breaks=seq(0,3.1, by=0.1))
	#get max and min

	densityH[i] <- max(histL[[i]]$density)
	histL[[i]]$densityScale <-histL[[i]]$density*(0.5/ densityH[i])
	maxS[i] <- round(max(sweRate$absRate[sweRate$gcID == i]),1)
}
#create box plot quantiles
quant <- list()
for(i in 1:5){
	quant[[i]] <- quantile(sweRate$absRate[sweRate$gcID == i], probs=c(0.025,0.25,0.50,0.75,0.975))
}

#break up names 
nameSplit1 <- character(0)
nameSplit2 <- character(0)
for(i in 1:5){
	nameSplit1[i] <- strsplit(intercept$name2[i], " ")[[1]][1]
	nameSplit2[i] <- strsplit(intercept$name2[i], " ")[[1]][2]
}
nameSplit2 <- ifelse(is.na(nameSplit2), " ", nameSplit2)

#plotting
xseq <- seq(1,10, by=2)
plotOrder <- c(1,4,2,3,5)

wd1 <- 18
hd1 <- 18
yl <- 0
yh <- 3.1

png(paste0(plotDI,"\\violin.png"), width = 20, height = 20, units = "cm", res=300)
	layout(matrix(c(1),ncol=1, byrow=TRUE), width=lcm(wd1),height=lcm(hd1))
	plot(c(0,1),c(0,1), xlim=c(0,12), ylim=c(yl,yh), axes=FALSE, type="n", xlab = " ", ylab= " ",
		xaxs="i", yaxs="i")
		for(j in 1:5){
		i <- plotOrder[j]
			polygon(c(xseq[j]+(0-histL[[i]]$densityScale[histL[[i]]$mids<=maxS[i]]), 
						rev(xseq[j]+(histL[[i]]$densityScale[histL[[i]]$mids<=maxS[i]]))),
					c(histL[[i]]$mids[histL[[i]]$mids<=maxS[i]],
						rev(histL[[i]]$mids[histL[[i]]$mids<=maxS[i]])), 
					lwd=0.75,  col=vegePallete3[i])
			arrows(	xseq[j],quant[[i]][1], xseq[j],quant[[i]][5], code=0, lwd=1)
			polygon(c(xseq[j]-0.15,xseq[j]-0.15,xseq[j]+0.15,xseq[j]+0.15),
				c(quant[[i]][2],quant[[i]][4],quant[[i]][4],quant[[i]][2]),
					border=NA, col=rgb(0.25,0.25,0.25,0.5))
					
			arrows(	xseq[j]-.15,quant[[i]][3], xseq[j]+.15,quant[[i]][3], code=0, lwd=4, col=vegePallete[i])		
		}	
		
	axis(1, xseq, rep(" ",length(xseq)),lwd.ticks=tlw)
	axis(2, seq(0,3.1, by=.2),rep(" ",length(seq(0,3.1, by=.2))),lwd.ticks=tlw)
	mtext(paste(nameSplit1[plotOrder]),at=xseq,side=1,line=1,cex=1)
	mtext(paste(nameSplit2[plotOrder]),at=xseq,side=1,line=2,cex=1)
	mtext(seq(0,3.1, by=.2), at=seq(0,3.1, by=.2), side=2, las=2, line=1, cex=1)
	mtext(expression(paste("Melt rate (cm day"^"-1",")")), side=2, line=3, cex=1.5)
	mtext("Landcover type", side=1, line=3, cex=1.5)
dev.off()




##################################################################################3
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 

swe.files <- list.files(pattern =".nc",path =paste0(swepath,"\\swe_mudryk_blended"),full.names=T)

#set up quantile plot
sweBreaks <- round(quantile(sweAll$swe,seq(0,1,length.out=11)),2)
sweBreaks[1] <- 0
# read one file in to use for reprojecting
pr <- brick(swe.files[29])
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')
worldmap <- map("world", ylim=c(40,90), fill=TRUE)


#world map
world <- project(matrix(c(worldmap$x,worldmap$y), ncol=2,byrow=FALSE),laea)


worldP <- Polygon(na.omit(matrix(c(world[,1],world[,2]), ncol=2, byrow=FALSE)))

worldSH <- readOGR("z:\\projects\\boreal_swe_depletion\\data\\world_cl.shp")
plot(worldSH)


worldr <- rasterize(worldSH, subset(pr,1))

plot(worldr)


pr3 <- mask(pr,worldr)

colP <- c(rgb(84,48,5, maxColorValue = 255),
				rgb(140,81,10, maxColorValue = 255),
				rgb(191,129,45, maxColorValue = 255),
				rgb(223,194,125, maxColorValue = 255),
				rgb(246,232,195, maxColorValue = 255),
				rgb(245,245,245, maxColorValue = 255),
				rgb(199,234,229, maxColorValue = 255),
				rgb(128,205,193, maxColorValue = 255),
				rgb(53,151,143, maxColorValue = 255),
				rgb(1,102,94, maxColorValue = 255),
				rgb(0,60,48, maxColorValue = 255))

land <- rgb(139,98,76,200, max=255)
				
colP2 <- c(land,
	rgb(255,247,251, maxColorValue = 255),
	rgb(236,231,242, maxColorValue = 255),
	rgb(208,209,230, maxColorValue = 255),
	rgb(166,189,219, maxColorValue = 255),
	rgb(116,169,207, maxColorValue = 255),
	rgb(54,144,192, maxColorValue = 255),
	rgb(5,112,176, maxColorValue = 255),
	rgb(4,90,141, maxColorValue = 255),
	rgb(2,56,88	, maxColorValue = 255))			
sweBreaks2 <- round(quantile(sweAll$swe,seq(0,1,length.out=11)),2)
sweBreaks2[1] <- 0

for(i in 1:152){
	png(paste0(plotDI,"\\swe_present\\swe00",i,".png"), width = 10, height = 10, units = "cm", res=300)

	#set up empty plot
	### plot 1 vegetation type ###
	par(mai=c(0.25,0.25,0.25,0.25))
	
	plot(subset(pr3,i), breaks=sweBreaks2,col=colP2,
			axes=FALSE, box=FALSE, legend=FALSE) 	
			
	legend("bottomright",paste("Day of year",i), bty="n")
	dev.off()
}
png(paste0(plotDI,"\\swe_legend.png"), width = 10, height = 10, units = "cm", res=300)
	par(mai=c(0.25,0.25,0.25,1))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(seq(0,1,length.out=11))-1)){
		polygon(c(0,0,1,1), 
			c(seq(0,1,length.out=11)[i],
			seq(0,1,length.out=11)[i+1],
			seq(0,1,length.out=11)[i+1],
			seq(0,1,length.out=11)[i]),
			col=colP2[i],border=NA)
	}
	axis(4,seq(0,1,length.out=11),sweBreaks2,cex.axis=2,las=2)	

dev.off()


#########################################################
#make plot of swe decline
#4004 is good
cellPlot <- 4004
yearPlot <- 2009


par(mai=c(1.5,1.5,1.5,1.5))
plot(sweAll$jday[sweAll$cell==cellPlot & sweAll$year == yearPlot], 
	sweAll$swe[sweAll$cell==cellPlot & sweAll$year == yearPlot],
	xlab=" ", ylab=" ",type="n",	axes=FALSE)
	polygon(c(130,130,180,180),
			c(0,.25,.25,0),
			col=rgb(205/255,79/255,57/255, .5),border=NA)
	points(sweAll$jday[sweAll$cell==cellPlot & sweAll$year == yearPlot], 
	sweAll$swe[sweAll$cell==cellPlot & sweAll$year == yearPlot],
	 pch=19)
axis(1,seq(30,180,by=30), cex.axis=1.5)
axis(2, seq(0,.25,by=0.05), las=2, cex.axis=1.5	)
mtext("Snow water equivalent (m)", side=2, line=5, cex=2)
mtext("Day of year", side=1, line=5, cex=2)	

#############################
#run on plot
plot(c(0,1),c(0,1))
legend("bottomright", paste(IDSglc$name2[c(1,2,4,3,5)]), col=vegePallete[c(1,2,4,3,5)],cex=2, lwd=2,lty=1, bty="n", bg="white")