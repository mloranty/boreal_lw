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
###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

modDir <- "z:\\projects\\boreal_swe_depletion\\analysis\\run6"

###############################################
### set up a dataframe with all of the      ###
### data and parameters by cell             ###
###############################################
sweCell <- unique(data.frame(cell=datSwe$cell,gcID=datSwe$gcID,pixID=datSwe$pixID,year=datSwe$year,
								y=datSwe$y.coord,x=datSwe$x.coord,
								vcf=datSwe$vcf,zone=datSwe$zone,dayMax=datSwe$dayMax))
								
								
#calculate the melt period		
#join midpoint into swe cell

sweCell2 <- join(sweCell,midOut,by=c("pixID","gcID","year"),type="left")

sweCell2$Mlength <- (round(sweCell2$Mean)-sweCell2$dayMax)*2
#calculate the end day
sweCell2$dayEnd <- sweCell2$dayMax+sweCell2$Mlength


#calculate average temperature to only be within the melt period
#and calculate the temperature in the three weeks before the melt period
meltTemp <- numeric(0)
onsetTemp <- numeric(0)


for(i in 1:dim(sweCell2)[1]){
	meltTemp[i] <- mean(sweAll$t.air[sweAll$pixID==sweCell2$pixID[i]&sweAll$gcID==sweCell2$gcID[i]&
						sweAll$year==sweCell2$year[i]&sweAll$jday>=sweCell2$dayMax[i]&sweAll$jday<=sweCell2$dayEnd[i]]-273.15)
	onsetTemp[i] <- mean(sweAll$t.air[sweAll$pixID==sweCell2$pixID[i]&sweAll$gcID==sweCell2$gcID[i]&
						sweAll$year==sweCell2$year[i]&sweAll$jday>=(sweCell2$dayMax[i]-7)&sweAll$jday<=sweCell2$dayMax[i]]-273.15)								
}

#create a data frame  to combine back into slope output
sweCell3 <- data.frame(pixID=sweCell2$pixID,cell=sweCell2$cell,gcID=sweCell2$gcID,year=sweCell2$year,
							vcf=sweCell2$vcf,zone=sweCell2$zone,dayMax=sweCell2$dayMax,Mlength=sweCell2$Mlength,
							dayEnd=sweCell2$dayEnd,meltTemp=meltTemp,onsetTemp=onsetTemp,
							x=sweCell2$x,y=sweCell2$y)

#now join to each output
b0All <- join(b0Out,sweCell3,by=c("year","gcID","pixID"),type="inner")


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


#jags regression
datalist <- list(Nobs= dim(b0All3)[1],
					maxD=b0All3$dayMax,
					b0=b0All3$Mean,
					glcIDM=b0All3$gcID,
					glcIDB=b0All3$gcID,
					TempAB=b0All3$meltTemp,
					CanopyB=b0All3$vcf,
					sweMaxB=b0All3$sweMax,
					TempAM=b0All3$onsetTemp,
					CanopyM=b0All3$vcf,
					Lat=b0All3$Lat,
					sig.modB=b0All3$SD,
					Nglc=dim(IDSglc)[1])

inits <- list(list(sig.vM=2,sig.vB=2),
				list(sig.vM=10,sig.vB=10),
				list(sig.vM=5,sig.vB=5))
				
parms <- c("betaM0","betaM1","betaM2","betaM3",
			"betaB0","betaB1","betaB2","betaB3",
			"mu.betaM0","mu.betaM1","mu.betaM2","mu.betaM3",
			"mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3",
			"sig.M0","sig.M1","sig.M2","sig.M3",
			"sig.B0","sig.B1","sig.B2","sig.B3",
			"sig.vB","sig.vM", "rep.mid","rep.b0")
			
	
curve.mod <- jags.model(file="c:\\Users\\hkropp\\Documents\\GitHub\\boreal_lw\\swe_model\\swe_curve_empirical_regression.r",
						data=datalist,n.adapt=5000,n.chains=3,inits=inits)
						
curve.sample <- coda.samples(curve.mod,variable.names=parms,n.iter=90000,thin=30)						
			
mcmcplot(curve.sample, parms=c("betaM0","betaM1","betaM2","betaM3",
			"betaB0","betaB1","betaB2","betaB3",
			"mu.betaM0","mu.betaM1","mu.betaM2","mu.betaM3",
			"mu.betaB0","mu.betaB1","mu.betaB2","mu.betaB3",
			"sig.M0","sig.M1","sig.M2","sig.M3",
			"sig.B0","sig.B1","sig.B2","sig.B3",
			"sig.vB","sig.vM"),dir=paste0(modDir,"\\history"))		


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
			
			
			
			
			