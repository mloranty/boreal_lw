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
library(rjags)
library(coda)
library(mcmcplots)
library(loo)
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

modDir <- "z:\\projects\\boreal_swe_depletion\\analysis\\run11"

###############################################
### set up a dataframe with all of the      ###
### data and parameters by cell             ###
###############################################
sweCell <- unique(data.frame(cell=datSwe$cell,gcID=datSwe$gcID,pixID=datSwe$pixID,year=datSwe$year,
								y=datSwe$y.coord,x=datSwe$x.coord,
								vcf=datSwe$vcf,zone=datSwe$zone,dayMax=datSwe$dayMax, newpixID=datSwe$newpixID))
								
								
#calculate the melt period		
#join midpoint into swe cell
sweCell2 <- join(sweCell,midOut,by=c("newpixID","gcID","year"),type="left")
colnames(halfOut)[1:7] <- paste0(colnames(halfOut)[1:7],"H")
sweCell3 <- join(sweCell2,halfOut,by=c("newpixID","gcID","year"),type="left")

#calculate the end day
sweCell3$dayEnd <- ifelse(round(sweCell3$Mean)+round(sweCell3$MeanH) > 182,
					182,round(sweCell3$Mean)+round(sweCell3$MeanH))

#### data filter #####
#there are some very fast melt periods
#where the standard deviation of the midpoint
#is within the range of the onset. Looking at the
#melt period won't be reliable. 					
sweCell3 <- sweCell3[round(sweCell3$Mean-sweCell3$SD)>sweCell3$dayMax,]

#calculate average temperature to only be within the melt period
#and calculate the temperature in the  week before the melt period

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
###############################################
### map results                             ###
###############################################


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
#plot each years curve on the map

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


#onset of melt
m2000 <- setValues(swe,bSwe[[1]]$dayMax)
m2001 <- setValues(swe,bSwe[[2]]$dayMax)
m2002 <- setValues(swe,bSwe[[3]]$dayMax)
m2003 <- setValues(swe,bSwe[[4]]$dayMax)
m2004 <- setValues(swe,bSwe[[5]]$dayMax)
m2005 <- setValues(swe,bSwe[[6]]$dayMax)
m2006 <- setValues(swe,bSwe[[7]]$dayMax)
m2007 <- setValues(swe,bSwe[[8]]$dayMax)
m2008 <- setValues(swe,bSwe[[9]]$dayMax)
m2009 <- setValues(swe,bSwe[[10]]$dayMax)
#plot each years curve on the map

plot(m2000)
plot(m2001)
plot(m2002)
plot(m2003)
plot(m2004)
plot(m2005)
plot(m2006)
plot(m2007)
plot(m2008)
plot(m2009)

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

#index for spatial random effect
cellDF <- unique(data.frame(cell=b0All4$cell, y=b0All4$y,x=b0All4$x,gcID=b0All4$gcID))
cellDF$cellID <- seq(1,nrow(cellDF))

#join back into b0All
b0All5 <- join(b0All4, cellDF, by=c("cell","x","y","gcID"), type="left")


#jags regression
datalist <- list(Nobs= dim(b0All5)[1],
					b0=b0All5$Mean,
					glcIDB=b0All5$gcID,
					TempAB=b0All5$meltTemp,
					CanopyB=b0All5$vcf,
					sweDay=b0All5$dayMax,
					GCyearB=b0All5$gcyearID,
					sig.modB=b0All5$SD,
					Ngcyear=dim(epsTable)[1],
					ygcIDB=epsTable$gcID,
					startb=startID,
					endb=endID,
					Nglc=dim(IDSglc)[1],
					cellID=b0All5$cellID,
					x=b0All5$x,
					y=b0All5$y,
					Ncell=nrow(cellDF))

inits <- list(list(sig.vB=2,sig.eb=rep(.5,dim(gcIndT)[1]),t.S=1,rhoS=.5,sigS=3,
					eps.s=rnorm(nrow(cellDF),0,.5)),
				list(sig.vB=10,sig.eb=rep(.6,dim(gcIndT)[1]),t.S=2,rhoS=.9,sigS=5,
						eps.s=rnorm(nrow(cellDF),-0.001,.1)),
				list(sig.vB=5,sig.eb=rep(.25,dim(gcIndT)[1]),t.S=3,rhoS=.6,sigS=1,
				eps.s=rnorm(nrow(cellDF),0.001,.25)))
				
parms <- c("betaB0S","betaB1","betaB2","betaB3",
			"mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3",
			"sig.B0","sig.B1","sig.B2","sig.B3",
			"sig.vB","rep.b0","eps.bS","sig.eb","Dsum","loglike","eps.sS","sigS","rhoS")
			
	
curve.mod <- jags.model(file="c:\\Users\\hkropp\\Documents\\GitHub\\boreal_lw\\swe_model\\swe_curve_empirical_regression.r",
						data=datalist,n.adapt=10000,n.chains=3,inits=inits)
						
curve.sample <- coda.samples(curve.mod,variable.names=parms,n.iter=50000,thin=25)						
			
mcmcplot(curve.sample, parms=c(
			"betaB0S","betaB1","betaB2","betaB3",
			"mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3",
			"sig.B0","sig.B1","sig.B2","sig.B3",
			"sig.vB","eps.bS","sig.eb","sigS","rhoS"),dir=paste0(modDir,"\\history"))		


#model output							   
mod.out <- summary(curve.sample)


write.table(mod.out$statistics,paste0(modDir,"\\curve_mod_stats.csv"),
			sep=",",row.names=TRUE)
write.table(mod.out$quantiles,paste0(modDir,"\\curve_mod_quant.csv"),
			sep=",",row.names=TRUE)	
			

#coda output
chain1<-as.matrix(curve.sample[[1]])
write.table(chain1,paste0(modDir,"\\chain1_coda.csv"), sep=",")
chain2<-as.matrix(curve.sample[[2]])
write.table(chain2,paste0(modDir,"\\chain2_coda.csv"), sep=",")
chain3<-as.matrix(curve.sample[[3]])
write.table(chain3,paste0(modDir,"\\chain3_coda.csv"), sep=",")					
			

###############################################
### read in regression results              ###
### check goodness of fit
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
#pull out slope rep
bRep <- datC[datC$parm=="rep.b0",]			
			
plot(b0All5$Mean,bRep$Mean)	
fit <- lm(bRep$Mean~	b0All5$Mean)	
summary(fit)			
abline(0,1,col="red",lwd=2) 
chains <- rbind(chain1,chain2,chain3)
llall <- chains[,gsub(dexps,"",colnames(chains))=="loglike"]
waic(llall)

datC[datC$parm=="Dsum",]

b0All5$residual <- b0All5$Mean-bRep$Mean

qqnorm(b0All5$residual)
qqline(b0All5$residual)


#check spatial patterns of residuals
#pull out by year
years <- seq(2000,2009)
b0All5L <- list()
for(i in 1:length(years)){
	b0All5L [[i]] <- join(sweCellDF,b0All5[b0All5$year==years[i],], by="cell",type="left")

}


#set into raster
rresid2000 <- setValues(swe,b0All5L [[1]]$residual)
rresid2001 <- setValues(swe,b0All5L [[2]]$residual)
rresid2002 <- setValues(swe,b0All5L [[3]]$residual)
rresid2003 <- setValues(swe,b0All5L [[4]]$residual)
rresid2004 <- setValues(swe,b0All5L [[5]]$residual)
rresid2005 <- setValues(swe,b0All5L [[6]]$residual)
rresid2006 <- setValues(swe,b0All5L [[7]]$residual)
rresid2007 <- setValues(swe,b0All5L [[8]]$residual)
rresid2008 <- setValues(swe,b0All5L [[9]]$residual)
rresid2009 <- setValues(swe,b0All5L [[10]]$residual)

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
