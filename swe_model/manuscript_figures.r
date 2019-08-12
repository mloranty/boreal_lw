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
colnames(midOut)[1:7] <-  paste0(colnames(midOut)[1:7],"M")
colnames(halfOut)[1:7] <- paste0(colnames(halfOut)[1:7],"H")
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

sweCell2 <- join(sweCell,midOut,by=c("newpixID","gcID","year"),type="left")
sweCell3 <- join(sweCell2,halfOut,by=c("newpixID","gcID","year"),type="left")

#### data filter #####
#there are some very fast melt periods
#where the standard deviation of the midpoint
#is within the range of the onset. Looking at the
#melt period won't be reliable. 					
sweCell3 <- sweCell3[round(sweCell3$MeanM-sweCell3$SDM)>sweCell3$dayMax,]

#parm all
parmAll <- join(b0Out,sweCell3,by=c("year","gcID","newpixID"),type="inner")

###############################################
### organize model data                     ###
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
betaMax <- datC[datC$parm=="betaB3",] 
betaTemp <- datC[datC$parm=="betaB1",] 
intC <- datC[datC$parm=="betaB0S",] 
#pull out slope rep
bRep <- datC[datC$parm=="rep.b0",]	


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
namesI <- c("Deciduous needleleaf boreal","Deciduous shrub tundra","Herbaceous tundra", "Evergreen needleleaf  boreal","Mixed leaf boreal" )


treePallete <- c(rgb(229,245,224,max=255),
				rgb(199,233,192,max=255),
				rgb(161,217,155,max=255),
				rgb(116,196,118,max=255),
				rgb(65,171,93,max=255),
				rgb(35,139,69,max=255),
				rgb(0,109,44,max=255),
				rgb(0,68,27,max=255))
vegePallete <- c(rgb(170/255,190/255,140/255),	
				rgb(60/255,60/255,110/255),
				rgb(130/255,160/255,190/255),
				rgb(50/255,80/255,10/255),
				rgb(250/255,120/255,80/255))
				
swePallete <- c(rgb(247,252,253,max=255),
				rgb(224,236,244,max=255),
				rgb(158,188,218,max=255),
				rgb(140,150,198,max=255),
				rgb(140,107,177,max=255),
				rgb(136,65,157,max=255),
				rgb(129,15,124,max=255),
				rgb(77,0,75,max=255))

hd <- 10
wd1 <- 10
wd2 <- 4
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
cxa <- 1.25

png(paste0(plotDI,"\\data_maps_vege.png"), width = 12, height = 5, units = "in", res=300)
	layout(matrix(seq(1,2),ncol=2), width=c(lcm(wd1),lcm(wd2)),height=lcm(hd))
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
	image(rasterVege,breaks=vegeBr,col=vegePallete,add=TRUE )
	mtext("A",at=4100000,side=2,line=pll, las=2,cex=mx)
	
	### plot 1 legend ###
	par(mai=c(0,.25,0,1))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(vegeBr)-1)){
		polygon(c(0,0,1,1), c(vegeBr[i]/(length(vegeBr)-1),vegeBr[i+1]/(length(vegeBr)-1),vegeBr[i+1]/(length(vegeBr)-1),vegeBr[i]/(length(vegeBr)-1)),col=vegePallete[i],border=NA)
	}
	axis(4,(vegeBr[1:5]/5)+.1,namesI,cex.axis=cxa,las=2)
	
dev.off()

png(paste0(plotDI,"\\data_maps_canopy.png"), width = 12, height = 5, units = "in", res=300)
	layout(matrix(seq(1,2),ncol=2), width=c(lcm(wd1),lcm(wd2)),height=lcm(hd))	
	
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
	image(rasterCanopy, breaks=canopyBr, col=treePallete, add=TRUE)
	mtext("B",at=4100000,side=2,line=pll, las=2,cex=mx)
	### plot 2 legend ###
	par(mai=c(0,.25,0,1))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(canopyBr)-1)){
		polygon(c(0,0,1,1), 
			c(canopyBr[i]/canopyBr[length(canopyBr)],canopyBr[i+1]/canopyBr[length(canopyBr)],canopyBr[i+1]/canopyBr[length(canopyBr)],canopyBr[i]/canopyBr[length(canopyBr)]),
			col=treePallete[i],border=NA)
	}
	axis(4,canopyBr/canopyBr[length(canopyBr)],canopyBr,cex.axis=cxa,las=2)	
	
dev.off()

png(paste0(plotDI,"\\data_maps_maxSwe.png"), width = 12, height = 5, units = "in", res=300)
	layout(matrix(seq(1,2),ncol=2), width=c(lcm(wd1),lcm(wd2)),height=lcm(hd))	

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
	image(rasterMaxMean,breaks=sweBr, col=swePallete, add=TRUE)
	mtext("C",at=4100000,side=2,line=pll, las=2,cex=mx)
	### plot 3 legend ###
	par(mai=c(0,.25,0,1))
	plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
	for(i in 1:(length(sweBr)-1)){
		polygon(c(0,0,1,1), 
			c(sweBr[i]/sweBr[length(sweBr)],sweBr[i+1]/sweBr[length(sweBr)],sweBr[i+1]/sweBr[length(sweBr)],sweBr[i]/sweBr[length(sweBr)]),
			col=swePallete[i],border=NA)
	}
	axis(4,sweBr/sweBr[length(sweBr)],sweBr,cex.axis=cxa,las=2)	
	
dev.off()


################################################################################
################################################################################
############### Figure 2. Plot regression                        ############### 
################################################################################
################################################################################


###############################################
### organize regression data                ###
###############################################

sweCell3$dayEnd <- ifelse(round(sweCell3$MeanM)+round(sweCell3$MeanH) > 182,
					182,round(sweCell3$MeanM)+round(sweCell3$MeanH))

#subset swe cell3 info
daysToJoin <- data.frame(pixID=sweCell3$pixID,year=sweCell3$year,gcID=sweCell3$gcID,
						dayMax=sweCell3$dayMax,dayEnd=sweCell3$dayEnd)

						
sweDaysJoin <- join(sweAll,daysToJoin, by=c("pixID","gcID","year"),type="inner")

sweMeltSub <- sweDaysJoin[sweDaysJoin$jday>=sweDaysJoin$dayMax&sweDaysJoin$jday<=sweDaysJoin$dayEnd,]
sweOnsetSub <- sweDaysJoin[sweDaysJoin$jday>=(sweDaysJoin$dayMax-7)&sweDaysJoin$jday<=sweDaysJoin$dayMax,]		

#now aggregate temperature
meltTempDF <- aggregate(sweMeltSub$t.air-273.15, by=list(sweMeltSub$pixID,sweMeltSub$gcID,sweMeltSub$year),FUN="mean",na.rm=TRUE)
colnames(meltTempDF) <- c("pixID","gcID","year","meltTemp")
onsetTempDF <- aggregate(sweOnsetSub$t.air-273.15, by=list(sweOnsetSub$pixID,sweOnsetSub$gcID,sweOnsetSub$year),FUN="mean",na.rm=TRUE)
colnames(onsetTempDF) <- c("pixID","gcID","year","onsetTemp")
#join back into sweCell3
sweCell3a <- join(sweCell3, meltTempDF, by=c("pixID","gcID","year"),type="left")	
sweCell3b <- join(sweCell3a, onsetTempDF, by=c("pixID","gcID","year"),type="left")	

#create a data frame  to combine back into slope output
sweCell4 <- data.frame(pixID=sweCell3b$pixID,newpixID=sweCell3b$newpixID,cell=sweCell3b$cell,gcID=sweCell3b$gcID,year=sweCell3b$year,
							vcf=sweCell3b$vcf,zone=sweCell3b$zone,dayMax=sweCell3b$dayMax,
							dayEnd=sweCell3b$dayEnd,meltTemp=sweCell3b$meltTemp,onsetTemp=sweCell3b$onsetTemp,
							x=sweCell3b$x,y=sweCell3b$y)

#now join to each output
b0All <- join(b0Out,sweCell4,by=c("year","gcID","newpixID"),type="inner")


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

for(i in 1:dim(yearDF)[1]){
	bSwe[[i]] <- join(sweCellDF,bOutL[[i]], by="cell",type="left")

}
#get lat long for each cell
#create a spatial points
sweSP <- SpatialPoints(unique(data.frame(x=b0All$x,y=b0All$y,cell=b0All$cell)), CRS(laea))
#transform for wgs lat long
sweSPr <- spTransform(sweSP, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

sweLL <- data.frame(sweSPr@coords)
colnames(sweLL) <- c("Lon","Lat","cell")
#join back into dataframe of results
b0All2 <- join(b0All,sweLL, by="cell",type="left")


#get unique swe max
swemax <- unique(data.frame(gcID=datSwe$gcID,cell=datSwe$cell,year=datSwe$year,sweMax=datSwe$sweMax))

#join

b0All3 <- join(b0All2,swemax,by=c("gcID","cell","year"),type="left")


#need to organize table for eps ids
epsTable <- unique(data.frame(gcID=b0All3$gcID,year=b0All3$year))
#this order will be by GCID
epsTable$gcyearID <- seq(1,dim(epsTable)[1])
#join back into b0
b0All4 <- join(b0All3,epsTable, by=c("gcID","year"),type="left")

#create index for averaging eps
gcIndT <- unique(data.frame(gcID=epsTable$gcID))
startID <- numeric(0)
endID <- numeric(0)

for(i in 1:dim(gcIndT)[1]){
		startID[i] <- head(which(epsTable$gcID==gcIndT$gcID[i]))[1]
		endID [i] <- tail(which(epsTable$gcID==gcIndT$gcID[i]))[6]
}




plot(b0All4$Mean,bRep$Mean, xlim=c(0,1) , ylim=c(0,1) )
abline(0,1,col="red")

b0All4$residual <- b0All4$Mean-bRep$Mean

qqnorm(b0All4$residual)
qqline(b0All4$residual)


#check spatial patterns of residuals
#pull out by year
years <- seq(2000,2009)
b0All4L <- list()
for(i in 1:length(years)){
	b0All4L [[i]] <- join(sweCellDF,b0All4[b0All4$year==years[i],], by="cell",type="left")

}


#set into raster
rresid2000 <- setValues(swe,b0All4L [[1]]$residual)
rresid2001 <- setValues(swe,b0All4L [[2]]$residual)
rresid2002 <- setValues(swe,b0All4L [[3]]$residual)
rresid2003 <- setValues(swe,b0All4L [[4]]$residual)
rresid2004 <- setValues(swe,b0All4L [[5]]$residual)
rresid2005 <- setValues(swe,b0All4L [[6]]$residual)
rresid2006 <- setValues(swe,b0All4L [[7]]$residual)
rresid2007 <- setValues(swe,b0All4L [[8]]$residual)
rresid2008 <- setValues(swe,b0All4L [[9]]$residual)
rresid2009 <- setValues(swe,b0All4L [[10]]$residual)

par(mfrow=c(2,5))
plot(rresid2000)
plot(rresid2001)
plot(rresid2002)
plot(rresid2003)
plot(rresid2004)
plot(rresid2005)
plot(rresid2006)
plot(rresid2007)
plot(rresid2008)
plot(rresid2009)

par(mfrow=c(3,2))
plot(b0All4$sweMax[b0All4$gcID==1],b0All4$Mean[b0All4$gcID==1],pch=19)
abline(intC$Mean[1],betaMax$Mean[1], col="red")
plot(b0All4$sweMax[b0All4$gcID==2],b0All4$Mean[b0All4$gcID==2],pch=19)
abline(intC$Mean[2],betaMax$Mean[2], col="red")
plot(b0All4$sweMax[b0All4$gcID==3],b0All4$Mean[b0All4$gcID==3],pch=19)
abline(intC$Mean[3],betaMax$Mean[3], col="red")
plot(b0All4$sweMax[b0All4$gcID==4],b0All4$Mean[b0All4$gcID==4],pch=19)
abline(intC$Mean[4],betaMax$Mean[4], col="red")
plot(b0All4$sweMax[b0All4$gcID==5],b0All4$Mean[b0All4$gcID==5],pch=19)
abline(intC$Mean[5],betaMax$Mean[5], col="red")
plot(seq(1,5), betaMax$Mean, xlab="glc", ylim=c(-0.8,0.1),pch=19)
arrows(seq(1,5), betaMax$X2.5.,seq(1,5), betaMax$X97.5.,code=0)

par(mfrow=c(3,2))
plot(b0All4$vcf[b0All4$gcID==1],b0All4$Mean[b0All4$gcID==1],pch=19)
abline(intC$Mean[1],betaCov$Mean[1], col="red")
plot(b0All4$vcf[b0All4$gcID==2],b0All4$Mean[b0All4$gcID==2],pch=19)
abline(intC$Mean[2],betaCov$Mean[2], col="red")
plot(b0All4$vcf[b0All4$gcID==3],b0All4$Mean[b0All4$gcID==3],pch=19)
abline(intC$Mean[3],betaCov$Mean[3], col="red")
plot(b0All4$vcf[b0All4$gcID==4],b0All4$Mean[b0All4$gcID==4],pch=19)
abline(intC$Mean[4],betaCov$Mean[4], col="red")
plot(b0All4$vcf[b0All4$gcID==5],b0All4$Mean[b0All4$gcID==5],pch=19)
abline(intC$Mean[5],betaCov$Mean[5], col="red")
plot(seq(1,5), betaCov$Mean, xlab="glc", ylim=c(-0.006,.003),pch=19)
arrows(seq(1,5), betaCov$X2.5.,seq(1,5), betaCov$X97.5.,code=0)

par(mfrow=c(3,2))
plot(b0All4$meltTemp[b0All4$gcID==1],b0All4$Mean[b0All4$gcID==1],pch=19)
abline(intC$Mean[1],betaTemp$Mean[1], col="red")
plot(b0All4$meltTemp[b0All4$gcID==2],b0All4$Mean[b0All4$gcID==2],pch=19)
abline(intC$Mean[2],betaTemp$Mean[2], col="red")
plot(b0All4$meltTemp[b0All4$gcID==3],b0All4$Mean[b0All4$gcID==3],pch=19)
abline(intC$Mean[3],betaTemp$Mean[3], col="red")
plot(b0All4$meltTemp[b0All4$gcID==4],b0All4$Mean[b0All4$gcID==4],pch=19)
abline(intC$Mean[4],betaTemp$Mean[4], col="red")
plot(b0All4$meltTemp[b0All4$gcID==5],b0All4$Mean[b0All4$gcID==5],pch=19)
abline(intC$Mean[5],betaTemp$Mean[5], col="red")
plot(seq(1,5), betaTemp$Mean, xlab="glc", ylim=c(0,.05),pch=19)
arrows(seq(1,5), betaTemp$X2.5.,seq(1,5), betaTemp$X97.5.,code=0)

plot(b0All4$meltTemp[b0All4$gcID==1],b0All4$vcf[b0All4$gcID==1])
plot(b0All4$meltTemp[b0All4$gcID==2],b0All4$vcf[b0All4$gcID==2])
plot(b0All4$meltTemp[b0All4$gcID==3],b0All4$vcf[b0All4$gcID==3])
plot(b0All4$meltTemp[b0All4$gcID==4],b0All4$vcf[b0All4$gcID==4])
plot(b0All4$meltTemp[b0All4$gcID==5],b0All4$vcf[b0All4$gcID==5])

plot(b0All4$meltTemp[b0All4$gcID==1],b0All4$sweMax[b0All4$gcID==1])
plot(b0All4$meltTemp[b0All4$gcID==2],b0All4$sweMax[b0All4$gcID==2])
plot(b0All4$meltTemp[b0All4$gcID==3],b0All4$sweMax[b0All4$gcID==3])
plot(b0All4$meltTemp[b0All4$gcID==4],b0All4$sweMax[b0All4$gcID==4])
plot(b0All4$meltTemp[b0All4$gcID==5],b0All4$sweMax[b0All4$gcID==5])

plot(b0All4$meltTemp[b0All4$gcID==1]*b0All4$vcf[b0All4$gcID==1],b0All4$Mean[b0All4$gcID==1],pch=19)
plot(b0All4$meltTemp[b0All4$gcID==2]*b0All4$vcf[b0All4$gcID==2],b0All4$Mean[b0All4$gcID==2],pch=19)
plot(b0All4$meltTemp[b0All4$gcID==3]*b0All4$vcf[b0All4$gcID==3],b0All4$Mean[b0All4$gcID==3],pch=19)
plot(b0All4$meltTemp[b0All4$gcID==4]*b0All4$vcf[b0All4$gcID==4],b0All4$Mean[b0All4$gcID==4],pch=19)
plot(b0All4$meltTemp[b0All4$gcID==5]*b0All4$vcf[b0All4$gcID==5],b0All4$Mean[b0All4$gcID==5],pch=19)

par(mfrow=c(3,2))
plot(b0All4$dayMax[b0All4$gcID==1],b0All4$Mean[b0All4$gcID==1])
plot(b0All4$dayMax[b0All4$gcID==2],b0All4$Mean[b0All4$gcID==2])
plot(b0All4$dayMax[b0All4$gcID==3],b0All4$Mean[b0All4$gcID==3])
plot(b0All4$dayMax[b0All4$gcID==4],b0All4$Mean[b0All4$gcID==4])
plot(b0All4$dayMax[b0All4$gcID==5],b0All4$Mean[b0All4$gcID==5])


plot(b0All4$dayMax[b0All4$gcID==1],b0All4$sweMax[b0All4$gcID==1])
plot(b0All4$dayMax[b0All4$gcID==2],b0All4$sweMax[b0All4$gcID==2])
plot(b0All4$dayMax[b0All4$gcID==3],b0All4$sweMax[b0All4$gcID==3])
plot(b0All4$dayMax[b0All4$gcID==4],b0All4$sweMax[b0All4$gcID==4])
plot(b0All4$dayMax[b0All4$gcID==5],b0All4$sweMax[b0All4$gcID==5])

plot(b0All4$dayMax[b0All4$gcID==1],b0All4$meltTemp[b0All4$gcID==1])
plot(b0All4$dayMax[b0All4$gcID==2],b0All4$meltTemp[b0All4$gcID==2])
plot(b0All4$dayMax[b0All4$gcID==3],b0All4$meltTemp[b0All4$gcID==3])
plot(b0All4$dayMax[b0All4$gcID==4],b0All4$meltTemp[b0All4$gcID==4])
plot(b0All4$dayMax[b0All4$gcID==5],b0All4$meltTemp[b0All4$gcID==5])


cor(b0All4$dayMax[b0All4$gcID==1],b0All4$meltTemp[b0All4$gcID==1])
cor(b0All4$dayMax[b0All4$gcID==2],b0All4$meltTemp[b0All4$gcID==2])
cor(b0All4$dayMax[b0All4$gcID==3],b0All4$meltTemp[b0All4$gcID==3])
cor(b0All4$dayMax[b0All4$gcID==4],b0All4$meltTemp[b0All4$gcID==4])
cor(b0All4$dayMax[b0All4$gcID==5],b0All4$meltTemp[b0All4$gcID==5])





par(mfrow=c(3,2))
plot(b0All4$dayMax[b0All4$gcID==1],b0All4$vcf[b0All4$gcID==1])
plot(b0All4$dayMax[b0All4$gcID==2],b0All4$vcf[b0All4$gcID==2])
plot(b0All4$dayMax[b0All4$gcID==3],b0All4$vcf[b0All4$gcID==3])
plot(b0All4$dayMax[b0All4$gcID==4],b0All4$vcf[b0All4$gcID==4])
plot(b0All4$dayMax[b0All4$gcID==5],b0All4$vcf[b0All4$gcID==5])


cor(b0All4$dayMax[b0All4$gcID==1],b0All4$vcf[b0All4$gcID==1])
cor(b0All4$dayMax[b0All4$gcID==2],b0All4$vcf[b0All4$gcID==2])
cor(b0All4$dayMax[b0All4$gcID==3],b0All4$vcf[b0All4$gcID==3])
cor(b0All4$dayMax[b0All4$gcID==4],b0All4$vcf[b0All4$gcID==4])
cor(b0All4$dayMax[b0All4$gcID==5],b0All4$vcf[b0All4$gcID==5])