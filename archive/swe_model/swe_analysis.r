##########################################################
#### Analysis of 9 years of swe depletion in boreal   ####
#### forests.                                         ####
##########################################################
#### inputs: daily swe data used in analysis: datSwe  ####
#### and all swe before subsetting melt period: allSwe####
#### table of melt period rate and other info: cellSwe####
#### zone ids: IDSglc                                 ####
##########################################################

###############################################
### read in swe depletion data              ###
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
library(rjags)
library(coda)
library(mcmcplots)
library(loo)
###############################################
### set up file paths                       ###
###############################################
swepath <- "C:\\Users\\hkropp\\Google Drive\\GIS"

modDir <- "C:\\Users\\hkropp\\Google Drive\\research\\projects\\boreal_swe\\boreal_swe_z\\analysis\\run22"

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


b2000 <- setValues(swe,meltSwe[[1]]$meltRateCM)
b2001 <- setValues(swe,meltSwe[[2]]$meltRateCM)
b2002 <- setValues(swe,meltSwe[[3]]$meltRateCM)
b2003 <- setValues(swe,meltSwe[[4]]$meltRateCM)
b2004 <- setValues(swe,meltSwe[[5]]$meltRateCM)
b2005 <- setValues(swe,meltSwe[[6]]$meltRateCM)
b2006 <- setValues(swe,meltSwe[[7]]$meltRateCM)
b2007 <- setValues(swe,meltSwe[[8]]$meltRateCM)
b2008 <- setValues(swe,meltSwe[[9]]$meltRateCM)
b2009 <- setValues(swe,meltSwe[[10]]$meltRateCM)
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
m2000 <- setValues(swe,meltSwe[[1]]$meltStart)
m2001 <- setValues(swe,meltSwe[[2]]$meltStart)
m2002 <- setValues(swe,meltSwe[[3]]$meltStart)
m2003 <- setValues(swe,meltSwe[[4]]$meltStart)
m2004 <- setValues(swe,meltSwe[[5]]$meltStart)
m2005 <- setValues(swe,meltSwe[[6]]$meltStart)
m2006 <- setValues(swe,meltSwe[[7]]$meltStart)
m2007 <- setValues(swe,meltSwe[[8]]$meltStart)
m2008 <- setValues(swe,meltSwe[[9]]$meltStart)
m2009 <- setValues(swe,meltSwe[[10]]$meltStart)
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
#set up data for plotting
tempMean <- seq(floor(range(cellSwe7$tair)[1]),ceiling(range(cellSwe7$tair)[2]), length.out=200)
CanopyMean <- seq(floor(range(cellSwe7$vcf)[1]),ceiling(range(cellSwe7$vcf)[2]), length.out=200)
SdayMean <- seq(floor(range(cellSwe7$meltStart)[1]),ceiling(range(cellSwe7$meltStart)[2]), length.out=200)
MaxMean <- seq(0,0.5, length.out=200)

#jags regression
datalist <- list(Nobs= dim(cellSwe7)[1],
					b0=cellSwe7$logAbsRate,
					glcIDB=cellSwe7$gcID,
					TempAB=cellSwe7$tair,
					CanopyB=cellSwe7$vcf,
					SweMax= cellSwe7$sweMax,
					sweDay=cellSwe7$meltStart,
					GCyearB=cellSwe7$gcyearID,
					Ngcyear=dim(epsTable)[1],
					ygcIDB=epsTable$gcID,
					startb=startID,
					endb=endID,
					Nglc=dim(IDSglc)[1],
					cellID=cellSwe7$cellID,
					Ncell=nrow(cellDF),
					TempMean=tempMean,
					CanopyMean=CanopyMean,
					SdayMean=SdayMean,
					MaxMean=MaxMean)


				
inits <- list(list(sig.eb=rep(.5,dim(gcIndT)[1]),sig.es=.2,
					eps.s=rnorm(nrow(cellDF),0,.5)),
				list(sig.eb=rep(.6,dim(gcIndT)[1]),sig.es=.5,
						eps.s=rnorm(nrow(cellDF),-0.001,.1)),
				list(sig.eb=rep(.25,dim(gcIndT)[1]),sig.es=1,
				eps.s=rnorm(nrow(cellDF),0.001,.25)))
				
parms <- c("betaB0S","betaB1","betaB2","betaB3","betaB4",
			"mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3","mu.betaB4",
			"sig.B0","sig.B1","sig.B2","sig.B3","sig.B4","trB0",
			"rep.b0","eps.bS","sig.eb","Dsum","loglike","eps.sS","sig.es",
			"mu.Temp","mu.Canopy","mu.Onset","mu.Max")
			
	
curve.mod <- jags.model(file="c:\\Users\\hkropp\\Documents\\GitHub\\boreal_lw\\swe_model\\swe_curve_empirical_regression.r",
						data=datalist,n.adapt=20000,n.chains=3,inits=inits)
						
curve.sample <- coda.samples(curve.mod,variable.names=parms,n.iter=600000,thin=300)						
			
mcmcplot(curve.sample, parms=c(
			"betaB0S","betaB1","betaB2","betaB3","betaB4",
			"mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3","mu.betaB4",
			"sig.B0","sig.B1","sig.B2","sig.B3","sig.B4",
			"sig.vB","eps.bS","sig.eb","sig.es"),dir=paste0(modDir,"\\history"))	

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
			
plot(cellSwe7$logAbsRate,bRep$Mean)	
fit <- lm(bRep$Mean~	cellSwe7$logAbsRate)	
summary(fit)			
abline(0,1,col="red",lwd=2) 
chains <- rbind(chain1,chain2,chain3)
llall <- chains[,gsub(dexps,"",colnames(chains))=="loglike"]
waic(llall)

datC[datC$parm=="Dsum",]

cellSwe7$residual <- cellSwe7$logAbsRate-bRep$Mean

qqnorm(cellSwe7$residual)
qqline(cellSwe7$residual)
plot(bRep$Mean,cellSwe7$residual)
par(mfrow=c(2,2))
plot(fit)


#check spatial patterns of residuals
#pull out by year
years <- seq(2000,2009)
cellSwe7L <- list()
for(i in 1:length(years)){
	cellSwe7L [[i]] <- join(sweCellDF,cellSwe7[cellSwe7$year==years[i],], by="cell",type="left")

}


#set into raster
rresid2000 <- setValues(swe,cellSwe7L [[1]]$residual)
rresid2001 <- setValues(swe,cellSwe7L [[2]]$residual)
rresid2002 <- setValues(swe,cellSwe7L [[3]]$residual)
rresid2003 <- setValues(swe,cellSwe7L [[4]]$residual)
rresid2004 <- setValues(swe,cellSwe7L [[5]]$residual)
rresid2005 <- setValues(swe,cellSwe7L [[6]]$residual)
rresid2006 <- setValues(swe,cellSwe7L [[7]]$residual)
rresid2007 <- setValues(swe,cellSwe7L [[8]]$residual)
rresid2008 <- setValues(swe,cellSwe7L [[9]]$residual)
rresid2009 <- setValues(swe,cellSwe7L [[10]]$residual)

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


library(geoR)
#check out semivariogram of residuals


bins <- 10
coords.mod <- as.matrix(cellSwe7[,c("x","y")]/1000)
max.dist <- 0.5 * max(dist(coords.mod))

v <- variog(coords = coords.mod, data = cellSwe7$residual, uvec = (seq(0,
 max.dist, length = bins)))


fit.v <- variofit(v, ini.cov.pars = c(0.015, 0.1),
  nugget =0.014, fix.nugget=FALSE,cov.model = "exponential",
 minimisation.function = "nls", weights = "equal")
plot(v)
lines(v)
summary(fit.v)




#pull out interaction beta
beta1 <- datC[datC$parm=="betaB1",] 
beta2 <- datC[datC$parm=="betaB2",] 
beta3 <- datC[datC$parm=="betaB3",] 
beta4 <- datC[datC$parm=="betaB4",] 
#intercept
betaIn <- datC[datC$parm=="betaB0S",] 
betaInE <- datC[datC$parm=="trB0",] 


