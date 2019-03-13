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
library(colorRamps)
library(viridis)
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

modDir <- "z:\\projects\\boreal_swe_depletion\\analysis\\run7"
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


#organize output by year
yearDF <- data.frame(year=unique(sweCell$year))
bOutL <- list()

for(i in 1:dim(yearDF)[1]){
	bOutL[[i]] <- b0All[b0All$year==yearDF$year[i],]

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


#get max swe
swemax <- unique(data.frame(gcID=datSwe$gcID,cell=datSwe$cell,year=datSwe$year,sweMax=datSwe$sweMax))

#pull out each year
swemaxL <- list()
	for(i in 1:dim(yearDF)[1]){
		swemaxL[[i]] <- join(sweCellDF,swemax[swemax$year==yearDF$year[i],], by="cell",type="left")
	}
#join into spatial data
max2000 <- setValues(swe,swemaxL[[1]]$sweMax)
max2001 <- setValues(swe,swemaxL[[2]]$sweMax)
max2002 <- setValues(swe,swemaxL[[3]]$sweMax)
max2003 <- setValues(swe,swemaxL[[4]]$sweMax)
max2004 <- setValues(swe,swemaxL[[5]]$sweMax)
max2005 <- setValues(swe,swemaxL[[6]]$sweMax)
max2006 <- setValues(swe,swemaxL[[7]]$sweMax)
max2007 <- setValues(swe,swemaxL[[8]]$sweMax)
max2008 <- setValues(swe,swemaxL[[9]]$sweMax)
max2009 <- setValues(swe,swemaxL[[10]]$sweMax)


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

worldmap <- map("world", ylim=c(0,90), fill=TRUE)

#world map
worldPR <- project(matrix(c(worldmap$x,worldmap$y), ncol=2,byrow=FALSE),laea)




#######################################
#####plot map  of slope of melt  ##### 
#######################################
tx <- 3
br1 <- 0
br2 <- 1
axisC <- 1.5
#set up breaks
breaks <- seq(br1,br2,by=.1)
breaks.lab <- seq(br1,br2,by=.2)
cols <- viridis_pal()(length(breaks)-1)


#2000
#set up empty plot
png(paste0(plotDir,"\\map_slope\\slope2000.png"),width = 480, height = 480, units = "px")
plot(b2000,col="white",breaks=breaks,  zlim=c(br1,br2), axes=FALSE,
		xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2000,col=cols,breaks=breaks, lab.breaks=breaks.lab, zlim=c(br1,br2),
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC),add=TRUE)		
mtext("2000", outer=TRUE,side=3,line=-3,cex=2)
dev.off()

#2001
png(paste0(plotDir,"\\map_slope\\slope2001.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2001,col="white",breaks=breaks, lab.breaks=breaks.lab, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
					axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2001,col=cols,breaks=breaks, zlim=c(br1,br2), lab.breaks=breaks,add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2001", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2002
png(paste0(plotDir,"\\map_slope\\slope2002.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2002,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2002,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2002", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2003
png(paste0(plotDir,"\\map_slope\\slope2003.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2003,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2003,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2003", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2004
png(paste0(plotDir,"\\map_slope\\slope2004.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2004,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2004,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2004", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2005
png(paste0(plotDir,"\\map_slope\\slope2005.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2005,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2005,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2005", outer=TRUE,side=3,line=-3,cex=tx)

dev.off()
#2006
png(paste0(plotDir,"\\map_slope\\slope2006.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2006,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2006,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2006", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2007
png(paste0(plotDir,"\\map_slope\\slope2007.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2007,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2007,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2007", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2008
png(paste0(plotDir,"\\map_slope\\slope2008.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2008,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2008,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2008", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2009
png(paste0(plotDir,"\\map_slope\\slope2009.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(b2009,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(b2009,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2009", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()



#plot all images
#read in all images

png(paste0(plotDir,"\\all_slope.png"),width=3000,height=6000)

	layout(matrix(seq(1,10),ncol=2,byrow=TRUE))
	for(i in 1:10){
		par(mai=c(0,0,0,0))
		plot(load.image(paste0(plotDir,"\\map_slope\\slope",i+1999,".png")),axes=FALSE)
	}	
	
		
dev.off()	

#######################################
#####plot3 map of day of swe max  ##### 
#######################################
tx <- 3
br1 <- 0
br2 <- .5
axisC <- 1.5
#set up breaks
breaks <- seq(br1,br2,by=.05)
breaks.lab <- seq(br1,br2,by=.1)
cols <-  viridis_pal()(length(breaks)-1)

#2000
#set up empty plot
png(paste0(plotDir,"\\map_max\\max2000.png"),width = 480, height = 480, units = "px")
plot(max2000,col="white",breaks=breaks,  zlim=c(br1,br2), axes=FALSE,
		xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2000,col=cols,breaks=breaks, lab.breaks=breaks.lab, zlim=c(br1,br2),
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC),add=TRUE)		
mtext("2000", outer=TRUE,side=3,line=-3,cex=2)
dev.off()

#2001
png(paste0(plotDir,"\\map_max\\max2001.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2001,col="white",breaks=breaks, lab.breaks=breaks.lab, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
					axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2001,col=cols,breaks=breaks, lab.breaks=breaks,add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2001", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()

#2002
png(paste0(plotDir,"\\map_max\\max2002.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2002,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2002,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2002", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2003
png(paste0(plotDir,"\\map_max\\max2003.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2003,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2003,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2003", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2004
png(paste0(plotDir,"\\map_max\\max2004.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2004,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2004,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2004", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2005
png(paste0(plotDir,"\\map_max\\max2005.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2005,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2005,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2005", outer=TRUE,side=3,line=-3,cex=tx)

dev.off()
#2006
png(paste0(plotDir,"\\map_max\\max2006.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2006,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2006,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2006", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2007
png(paste0(plotDir,"\\map_max\\max2007.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2007,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2007,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2007", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2008
png(paste0(plotDir,"\\map_max\\max2008.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2008,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2008,,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2008", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()
#2009
png(paste0(plotDir,"\\map_max\\max2009.png"),width = 480, height = 480, units = "px")
#set up empty plot
plot(max2009,col="white",breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),axes=FALSE,xlab=" ", ylab=" ",
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000),
			border=NA, col=rgb(180/255,205/255,205/255,.5))
#boundaries
points(worldPR, type="l", lwd=2, col="grey65")
#continent color
polygon(c(worldPR[,1],rev(worldPR[,1])), c(worldPR[,2],rev(worldPR[,2])),col="cornsilk2",border=NA)	
plot(max2009,col=cols,breaks=breaks, lab.breaks=breaks, zlim=c(br1,br2),add=TRUE,
		axis.args=list(at=breaks.lab,
						labels=breaks.lab,
						cex.axis=axisC))	
mtext("2009", outer=TRUE,side=3,line=-3,cex=tx)
dev.off()



#plot all images
#read in all images

png(paste0(plotDir,"\\all_max.png"),width=4000,height=5000)

	layout(matrix(seq(1,10),ncol=3,byrow=TRUE))
	for(i in 1:10){
		par(mai=c(0,0,0,0))
		plot(load.image(paste0(plotDir,"\\map_max\\max",i+1999,".png")),axes=FALSE)
	}	
	
	
		
dev.off()	



###############################################
### read in regression results              ###
###############################################

#read in model output
datS <- read.csv(paste0(modDir,"\\curve_mod_stats.csv"))
datQ <- read.csv(paste0(modDir,"\\curve_mod_quant.csv"))

#combine data frames
datC <- cbind(datS,datQ)
#pull out parameter names
dexps<-"\\[*[[:digit:]]*\\]"
datC$parm <- gsub(dexps,"",rownames(datC))

#pull out betaB2
betaCov <- datC[datC$parm=="betaB2",] 
betaCov$gcID <- seq(1,5)

#names for plotting
IDSglc$plotNames <- c("Needleleaf deciduous boreal","Deciduous shrub tundra","Herbaceous tundra","Evergreen needleaf boreal","Mixed boreal")

#join two together
betaCov2 <- join(betaCov,IDSglc,by="gcID",type="left")
betaCov2$plotOrder <- c(2,4,5,1,3)
betaCov3 <- betaCov2[order(betaCov2$plotOrder),]
wd <- 20
hd <- 20
xseq <- c(1,2.5,4,5.5,7)

###############################################
### plot melt slope with canopy cover       ###
###############################################


#make a plot of slope with tree cover
png(paste0(plotDir,"\\slope_tree_cover.png"),width = 800, height = 800, units = "px")
	#layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
	par(mai=c(5,3,1,1))
		plot(c(0,1),c(0,1), type="n", xlim=c(0,8),ylim=c(-.006,.001), xlab=" ", ylab=" ", xaxs="i",yaxs="i", axes=FALSE)
			
			points(c(0,8),c(0,0),type="l",lwd=3, lty=3,col="grey75")
			for(i in 1:5){
				polygon(c(xseq[i]-0.5,xseq[i]-0.5,xseq[i]+0.5,xseq[i]+0.5),
					c(betaCov3$X25.[i],betaCov3$X75.[i],betaCov3$X75.[i],betaCov3$X25.[i]),
						col=rgb(205/255,79/255,57/255,.85),border=NULL)
			}
			arrows(xseq-.5,betaCov3$Mean,xseq+.5,betaCov3$Mean,lwd=3,code=0)
			arrows(xseq,betaCov3$X2.5.,xseq,betaCov3$X97.5.,code=0,lwd=2)
			axis(2,seq(-.006,0.001,by=0.001), label=rep(" ",length(seq(-.006,0.001,by=0.001))), cex.axis=2)
			mtext(seq(-.006,0.001,by=0.001),at=seq(-.006,0.001,by=0.001),side=2, line=1, cex=2,las=2)
			axis(1,xseq, label=rep(" ",length(xseq)), cex.axis=2)
			text(xseq, rep(-.0065,length(xseq)),betaCov3$plotNames,srt=35,adj=1,cex=2,xpd=TRUE)
			mtext("Change in melt rate" , side=2,line=10,cex=3)
			mtext("per increase in % tree cover", side=2,line=7,cex=3)
			mtext("Landcover type", side=1,line=15,cex=3)
			
dev.off()		
###############################################
### plot melt slope means                   ###
###############################################
#make a plot of melt rate average
xseq <- rep(seq(1,7,by=1.5),times=10)+rep(seq(0,90,by=10),each=5)
IDSglc$plotOrder <- c(2,4,5,1,3)				
muB2 <- join(muB0Out, IDSglc, by="gcID",type="left")	
		
		

muB3 <-muB2[order(muB2$year,muB2$plotOrder),]	

plotNamesDF <-	data.frame(plotOrder=seq(1,5), plotNames=	
							c("Evergreen needleaf boreal","Needleleaf deciduous boreal","Mixed boreal","Deciduous shrub tundra","Herbaceous tundra"))

muB3 <- join(muB3,plotNamesDF,by="plotOrder",type="left")							
#indicate vegeclass colors
cols <- c(rgb(50/255,80/255,10/255),
rgb(170/255,190/255,140/255),
rgb(240/255,240/255,50/255),
rgb(250/255,230/255,140/255),
rgb(0/255,110/255,130/255))

colst <- c(rgb(50/255,80/255,10/255,.5),
rgb(170/255,190/255,140/255,.5),
rgb(240/255,240/255,50/255,.5),
rgb(250/255,230/255,140/255,.5),
rgb(0/255,110/255,130/255,.5))

muB3$cols <- rep(cols,times=10)		
#make a plot of slope with tree cover
png(paste0(plotDir,"\\glc_slope.png"),width = 1000, height = 800, units = "px")
	#layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
	par(mai=c(2,2,1,1))
	plot(c(0,1),c(0,1), type="n", xlim=c(0,98),ylim=c(0,0.6), xlab=" ", ylab=" ", xaxs="i",yaxs="i", axes=FALSE)
		points(xseq[muB3$plotOrder==1],muB3$Mean[muB3$plotOrder==1],type="l",lwd=2,col=colst[1])
		points(xseq[muB3$plotOrder==2],muB3$Mean[muB3$plotOrder==2],type="l",lwd=2,col=colst[2])
		points(xseq[muB3$plotOrder==3],muB3$Mean[muB3$plotOrder==3],type="l",lwd=2,col=colst[3])
		points(xseq[muB3$plotOrder==4],muB3$Mean[muB3$plotOrder==4],type="l",lwd=2,col=colst[4])
		points(xseq[muB3$plotOrder==5],muB3$Mean[muB3$plotOrder==5],type="l",lwd=2,col=colst[5])
		for(i in 1:50){	
		polygon(c(xseq[i]-0.5,xseq[i]-0.5,xseq[i]+0.5,xseq[i]+0.5),
					c(muB3$p25[i],muB3$p75[i],muB3$p75[i],muB3$p25[i]),
						col=muB3$cols[i],border=NA)
			}
		arrows(xseq-0.5,muB3$Mean,xseq+.5,muB3$Mean,code=0,lwd=2)	
		arrows(xseq,muB3$p2.5,xseq,muB3$p97.5,code=0,lwd=1.5)
		
		axis(2,seq(0,0.7,by=0.1),cex.axis=2,las=2)
		axis(1, xseq[seq(3,48,by=5)],lab=seq(2000,2009),cex.axis=2)
		legend(3,.61,paste(muB3$plotNames[1:5]),fill=muB3$cols[1:5],bty="n",cex=1.5)
		mtext("Mean melt rate", side=2, line=5, cex=2)
		mtext("Year", side=1, line=5, cex=2)
dev.off()	


###############################################
### show only two years                     ###
###############################################



xseq <- c(seq(1,7,by=1.5),(seq(1,7,by=1.5)+10))
			
		
muB4 <- muB3[muB3$year<=2002&muB3$year>=2001,]			
#make a plot of slope with tree cover
png(paste0(plotDir,"\\glc_slope_sub.png"),width = 1000, height = 800, units = "px")
	#layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
	par(mai=c(2,2,1,1))
	plot(c(0,1),c(0,1), type="n", xlim=c(0,19),ylim=c(0,0.62), xlab=" ", ylab=" ", xaxs="i",yaxs="i", axes=FALSE)
		for(i in 1:50){	
		polygon(c(xseq[i]-0.5,xseq[i]-0.5,xseq[i]+0.5,xseq[i]+0.5),
					c(0,muB4$Mean[i],muB4$Mean[i],0),
						col=muB4$cols[i],border=NA)
			}
			
		arrows(xseq,muB4$p2.5,xseq,muB4$p97.5,code=0,lwd=1.5)
		axis(2,seq(0,0.7,by=0.1),cex.axis=2.5,las=2)
		axis(1, c(-1,xseq[c(3,8)],20),lab=c(" ",seq(2001,2002)," "),cex.axis=2.5)
		legend("topright",paste(muB4$plotNames[1:5]),fill=muB4$cols[1:5],bty="n",cex=2)	
			mtext("Year", side=1,line=4,cex=3)
			mtext("Rate of snow melt", side=2,line=5,cex=3)
dev.off()	
	