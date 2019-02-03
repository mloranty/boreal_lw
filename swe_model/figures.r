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
library(sp)
library(rjags)
library(coda)
library(mcmcplots)
library(maps)
library(imager)
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

modDir <- "z:\\projects\\boreal_swe_depletion\\analysis\\run1"
plotDir <- "z:\\projects\\boreal_swe_depletion\\figures\\figures"

###############################################
### set up a dataframe with all of the      ###
### data and parameters by cell             ###
###############################################
sweCell <- unique(data.frame(cell=datSwe$cell,gcID=datSwe$gcID,pixID=datSwe$pixID,year=datSwe$year,
								y=datSwe$y.coord,x=datSwe$x.coord,
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


#get lat long for each cell
#create a spatial points
sweSP <- SpatialPoints(unique(data.frame(x=b0All$x,y=b0All$y,cell=b0All$cell)), CRS(laea))
#transform for wgs lat long
sweSPr <- spTransform(sweSP, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

sweLL <- data.frame(sweSPr@coords)
colnames(sweLL) <- c("Lon","Lat","cell")
#join back into dataframe of results
b0All2 <- join(b0All,sweLL, by="cell",type="left")
midAll2 <- join(midAll,sweLL, by="cell",type="left")

###############################################
### panel of results                        ###
###############################################
#get map data
world <- readOGR("c:\\Program Files (x86)\\ArcGIS\\Desktop10.5\\ArcGlobeData\\continent.shp")
worldP <- spTransform(world,laea)
plot(mid2000)

#######################################
#####world data                   ##### 
#######################################

worldmap <- map("world", ylim=c(50,90), fill=TRUE)

#world map
worldPR <- project(matrix(c(worldmap$x,worldmap$y), ncol=2,byrow=FALSE),laea)

#######################################
#####plot1 map of day of half melt##### 
#######################################
tx <- 3
br1 <- 30
br2 <- 170
#set up breaks
breaks <- seq(br1,br2,by=5)
breaks.lab <- seq(br1,br2,by=10)
cols <- heat.colors(breaks)

#2000
#set up empty plot
png(paste0(plotDir,"\\map_mid\\mid2000.png"),width = 480, height = 480, units = "px")
plot(mid2000,col="white",breaks=breaks, lab.breaks=breaks.lab, zlim=c(br1,br2), axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2000,col=cols,breaks=breaks, lab.breaks=breaks.lab, zlim=c(br1,br2),cex.lab=2,add=TRUE)		
mtext("2000", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()

#2001
png(paste0(plotDir,"\\map_mid\\mid2001.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2001,col="white",breaks=breaks, lab.breaks=breaks.lab, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2001,col=cols,breaks=breaks, lab.breaks=breaks,add=TRUE)	
mtext("2001", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2002
png(paste0(plotDir,"\\map_mid\\mid2002.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2002,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2002,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2002", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2003
png(paste0(plotDir,"\\map_mid\\mid2003.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2003,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2003,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2003", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2004
png(paste0(plotDir,"\\map_mid\\mid2004.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2004,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2004,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2004", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2005
png(paste0(plotDir,"\\map_mid\\mid2005.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2005,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2005,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2005", outer=TRUE,side=3,line=-3,cex=tx)

dev.off()
#2006
png(paste0(plotDir,"\\map_mid\\mid2006.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2006,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2006,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2006", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2007
png(paste0(plotDir,"\\map_mid\\mid2007.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2007,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2007,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2007", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2008
png(paste0(plotDir,"\\map_mid\\mid2008.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2008,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2008,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2008", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2009
png(paste0(plotDir,"\\map_mid\\mid2009.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(mid2009,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ")
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col="cornsilk2",border=NA)	
plot(mid2009,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE)	
mtext("2009", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()



#plot all images
#read in all images

png(paste0(plotDir,"\\all_mid.png"),width=6000,height=2000)

	layout(matrix(seq(1,10),ncol=5,byrow=TRUE))
	for(i in 1:10){
		par(mai=c(0,0,0,0))
		plot(load.image(paste0(plotDir,"\\map_mid\\mid",i+1999,".png")),axes=FALSE)
	}	
	
		
dev.off()	