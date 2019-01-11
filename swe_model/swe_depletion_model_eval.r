#######################################################
# Script for reading in a model run on a blended swe  #
# product to look at differences in timing of         #
# swe depletion                                       #
#######################################################
library(plyr)
library(coda)
library(mcmcplots)
#######################################################
# data info                                           #
#######################################################

#linux =1 or windows =2
runOS <- 2
#linux data directory first option, windows second optioon
DDdir <- c("/home/hkropp/boreal/data",
				"z:\\projects\\boreal_swe_depletion\\data")

#######################################################
# file name information                               #
#######################################################

#read in rep information


#create a vector with all of the parent directories
dirP <- c("z:\\projects\\boreal_swe_depletion\\model\\run11\\run1")

outD <- "z:\\projects\\boreal_swe_depletion\\model\\run11\\eval"
#get a list of all of the files


	dirA <- list.dirs(paste0(dirP), full.names=FALSE)
	#remove parent directory
	dirA <- dirA[-1]
	dimA <- length(dirA)
	datA <- data.frame(files=dirA,dirP=rep(dirP,dimA))



#pull out identifying info

for(i in 1:dim(datA)[1]){
		#vegetation land class
		datA$glc[i] <- as.numeric(strsplit(as.character(datA$files[i]),"_")[[1]][2])
		#year
		datA$year[i] <- as.numeric(strsplit(as.character(datA$files[i]),"_")[[1]][4])
		#chain
		datA$chain[i] <- as.numeric(strsplit(as.character(datA$files[i]),"_")[[1]][6])
}	
whichrep <- read.csv("z:\\projects\\boreal_swe_depletion\\data\\rep_subID.csv")
#turn chains into a list
chain1 <- datA[datA$chain==1,1:4]
colnames(chain1)[1:2] <- paste0(colnames(chain1)[1:2],"1") 

chain2 <- datA[datA$chain==2,1:4]
colnames(chain2)[1:2] <- paste0(colnames(chain2)[1:2],"2")

chain3 <- datA[datA$chain==3,1:4]
colnames(chain3)[1:2] <- paste0(colnames(chain3)[1:2],"3")

chains <- list(chain1,chain2,chain3)

chainDF <- join_all(chains,by=c("glc","year"),type="inner")



#get a list of all of the names of the variables
parms <- list.files(paste0(chainDF$dirP1[1],"\\",chainDF$files1[1]))

#read in coda
#for(i 1:dim(chainDF)[1]){
b0Out1 <- data.frame()
b0Out2 <- data.frame()
b0Out3 <- data.frame()
b0mcmc <- mcmc.list()
b0.diag <- list()

for(i in 1:dim(chainDF)[1]){

	#read in output
	b0Out1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\b0_out.csv"))
	b0Out2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\b0_out.csv"))
	b0Out3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\b0_out.csv"))
	colnames(b0Out1) <- paste0("b0_",seq(1,dim(b0Out1)[2]))
	colnames(b0Out2) <- paste0("b0_",seq(1,dim(b0Out2)[2]))
	colnames(b0Out3) <- paste0("b0_",seq(1,dim(b0Out3)[2]))
	#turn into mcmc
	b0Out1 <- as.mcmc(b0Out1)
	b0Out2 <- as.mcmc(b0Out2)
	b0Out3 <- as.mcmc(b0Out3)
	b0mcmc <- mcmc.list(b0Out1,b0Out2,b0Out3)
	#make plots
	dir.create(paste0(outD,"\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
	mcmcplot(b0mcmc,dir=paste0(outD,"\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
	b0.diag[[i]] <- gelman.diag(b0mcmc)
}
	

midOut1 <- data.frame()
midOut2 <- data.frame()
midOut3 <- data.frame()
midmcmc <- mcmc.list()
mid.diag <- list()
for(i in 1:dim(chainDF)[1]){
	#read in output
	midOut1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\mid0_out.csv"))
	midOut2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\mid0_out.csv"))
	midOut3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\mid0_out.csv"))
	colnames(midOut1) <- paste0("mid0_",seq(1,dim(midOut1)[2]))
	colnames(midOut2) <- paste0("mid0_",seq(1,dim(midOut2)[2]))
	colnames(midOut3) <- paste0("mid0_",seq(1,dim(midOut3)[2]))
	#turn into mcmc
	midOut1 <- as.mcmc(midOut1)
	midOut2 <- as.mcmc(midOut2)
	midOut3 <- as.mcmc(midOut3)
	midmcmc <- mcmc.list(midOut1,midOut2,midOut3)
	#make plots
	#dir.create(paste0(outD,"\\mid\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
	#mcmcplot(midmcmc,dir=paste0(outD,"\\mid\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
	mid.diag[[i]] <- gelman.diag(midmcmc)
}

#create a dataframe indicating the convergence of each gridcell and glc
midCheck <- list()
b0Check <- list()
for(i in 1:dim(chainDF)[1]){

		midCheck[[i]] <- data.frame(conv=ifelse(mid.diag[[i]]$psrf[,1]>1.2,1,0),
						gridID=seq(1,dim(mid.diag[[i]]$psrf)[1]),
						glc=rep(chainDF$glc[i],dim(mid.diag[[i]]$psrf)[1]),
						year=rep(chainDF$year[i],dim(mid.diag[[i]]$psrf)[1]))
						
		b0Check[[i]] <- data.frame(conv=ifelse(b0.diag[[i]]$psrf[,1]>1.2,1,0),
						gridID=seq(1,dim(b0.diag[[i]]$psrf)[1]),
						glc=rep(chainDF$glc[i],dim(b0.diag[[i]]$psrf)[1]),
						year=rep(chainDF$year[i],dim(b0.diag[[i]]$psrf)[1]))
}
#

#turn into a data frame
midConv <- ldply(midCheck,data.frame)
b0Conv <- ldply(b0Check,data.frame)

#######################################################
# read in and filter data                             #
#######################################################

#######################################################
# read in and filter data                             #
#######################################################
#read in data files
if(runOS==1){
	dat.swe <- read.csv(paste0(DDdir[1], "/swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[1], "/glc50_table.csv"))
	whichrep <- read.csv(paste0(DDdir[1],"/rep_subID.csv"))

}else{

	dat.swe <- read.csv(paste0(DDdir[2],"\\swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[2], "/glc50_table.csv"))
	whichrep <- read.csv(paste0(DDdir[2],"\\rep_subID.csv"))
}


print("finish reading in data")

##########################
##### Subset point 1 #####
##########################
#only focus on 2000-2009 for now
dat.swe <- dat.swe[dat.swe$year<=2009&dat.swe$year>=2000,]




#calculate proportion of land cover
dat.swe$glc1.p <- dat.swe$glc1f/3136
hist(dat.swe$glc1.p )
dat.swe$glc2.p <- dat.swe$glc2f/3136
hist(dat.swe$glc2.p )

#just use first glc class 				
dat.swe$zone <- dat.swe$glc1


##########################
##### Filter point 1 #####
##########################
#start with the simplest classification of glc 
#>50% is dominated by the landcover class
dat.swe1 <- dat.swe[dat.swe$glc1.p>=.5,]


##########################
##### Filter point 2 #####
##########################
#filter out land cover types not of interest
#don't include glc 3,9,10,14,15,16,18,19,20,21,22
dat.swe2 <- dat.swe1[dat.swe1$zone!=2&dat.swe1$zone!=11&dat.swe1$zone!=3&dat.swe1$zone!=15&dat.swe1$zone!=16&dat.swe1$zone!=10&dat.swe1$zone!=14&
				dat.swe1$zone!=9&dat.swe1$zone<18,]

##########################
##### Filter point 3 #####
##########################				
				
#filter out areas with no snow extent in spring
sweMax <- aggregate(dat.swe2$swe,by=list(dat.swe2$cell,dat.swe2$year),FUN="max")
colnames(sweMax) <- c("cell","year","sweMax")
sweMin <- aggregate(dat.swe2$swe,by=list(dat.swe2$cell,dat.swe2$year),FUN="min")
colnames(sweMin) <- c("cell","year","sweMin")
sweMax$Diff <- sweMax$sweMax-sweMin$sweMin
te <- hist(sweMin$sweMin, breaks=seq(0,0.5, by=0.01))
te2 <- hist(sweMax$Diff, breaks=seq(0,0.55, by=0.01))
te3 <- hist(sweMax$sweMax, breaks=seq(0,0.6, by=0.01))
#vast majority of minimum swe 97% is below 0.02 in a year
#exclude any sites and years where the maximum swe does not exceed 0.04 in the year
sweMaxF <- sweMax[sweMax$sweMax >=0.01,]
#join back to swe to filter
dat.swe3 <- join(dat.swe2, sweMaxF, by=c("cell","year"), type="inner")


#######################################################
# organize data for model run                         #
#######################################################
#create a table of identifiers
#unique glc

IDSglc <- unique(data.frame(zone=dat.swe3$zone))
IDSglc <- join(IDSglc, dat.glc[,1:2], by="zone", type="left")
IDSglc$gcID <- seq(1,dim(IDSglc)[1])


#join glc ID into dataframe
dat.swe4 <- join(dat.swe3,IDSglc, by="zone", type="left")

#normalize swe
#calculate percent
dat.swe4$sweP <- dat.swe4$swe/dat.swe4$sweMax

#round swe for 20% of peak to 1 and 20% of low to zero
dat.swe4$sweN <- ifelse(dat.swe4$sweP>=0.8,1,
					ifelse(dat.swe4$sweP<=0.2,0,(dat.swe4$sweP-0.2)/0.6))
							
					
#get unique pixel id in each glc for parameter id

pixID <- unique(data.frame(cell=dat.swe4$cell,year=dat.swe4$year, gcID=dat.swe4$gcID))

#want each cell in each year, gcID that will be the subset to run each model
#create year x gcID dataframe
gcYearID <- unique(data.frame(year=dat.swe4$year, gcID=dat.swe4$gcID))

#organize gcYear by the count of cells 
pixCountList <- list()
pixGLCTemp <- numeric(0)
for(i in 1:dim(gcYearID)[1]){
	pixCountList[[i]] <- pixID[pixID$gcID==gcYearID$gcID[i]&pixID$year==gcYearID$year[i],]
	pixGLCTemp[i] <- dim(pixCountList[[i]])[1]
}

gcYearID$cellCount <- pixGLCTemp

#organize gcYearID by cell count
gcYearID <- gcYearID[order(gcYearID$cellCount),]

gcYearID$gcYearID <- seq(1,dim(gcYearID)[1])

#subset into each glc xYear
pixList <- list()
pixGLC <- numeric(0)
for(i in 1:dim(gcYearID)[1]){
	pixList[[i]] <- pixID[pixID$gcID==gcYearID$gcID[i]&pixID$year==gcYearID$year[i],]
	pixList[[i]]$pixID <- seq(1,dim(pixList[[i]])[1])
	pixList[[i]]$gcYearID <- rep(gcYearID$gcYearID[i], dim(pixList[[i]])[1])
	pixGLC[i] <- dim(pixList[[i]])[1]
}

pixJ <- ldply(pixList,data.frame)

#join back into swe
dat.swe5 <- join(dat.swe4,pixJ, by=c("cell","year","gcID"), type="left")


print("finish data organize")
#######################################################
# subset swe to only use the swe                      #
# after swe max is reached                            #
#######################################################
#get the day that the final max occurs
maxTemp <- numeric(0)
maxN <- numeric(0)
#get the final swe max time
for(i in 1:dim(gcYearID)[1]){
	for(j in 1:dim(pixList[[i]])[1]){
		maxTemp <- which(dat.swe5$pixID==pixList[[i]]$pixID[j]&dat.swe5$year==pixList[[i]]$year[j]& dat.swe5$gcID==pixList[[i]]$gcID[j]&dat.swe5$sweN==1) 
		
		maxN[j] <- tail(maxTemp, n=1)
	}
	pixList[[i]]$finalMax <- maxN
}

pixJ2 <- ldply(pixList,data.frame)

pixJ2$dayMax <- dat.swe5$jday[pixJ2$finalMax]


dat.swe6 <- join(dat.swe5,pixJ2, by=c("cell","year","gcID","pixID","gcYearID"), type="left")

dat.swe7 <- dat.swe6[dat.swe6$jday>=dat.swe6$dayMax,]




#pull out which rows each gc is related
#sweRows <- list()
#sweDims <- numeric(0)
#for(i in 1:dim(gcYearID)[1]){
#	sweRows[[i]] <- which(dat.swe7$gcID==gcYearID$gcID[i]&dat.swe7$year==gcYearID$year[i])
#	sweDims[i] <- length(sweRows[[i]])
#}
#find out which 
#whichrep <- list()
#for(i in 1:dim(gcYearID)[1]){
#	if(sweDims[i]>5000){
#		whichrep[[i]] <- data.frame(gcID=rep(gcYearID$gcID[i], each=5000),year=rep(gcYearID$year[i], each=5000), rows=sample(sweRows[[i]],5000))
#	}else{
#		whichrep[[i]] <-data.frame(gcID=rep(gcYearID$gcID[i], each=sweDims[i]),year=rep(gcYearID$year[i], each=sweDims[i]), rows=sweRows[[i]])
#	}
#}
#whichrep <- ldply(whichrep,data.frame)
#write.table(whichrep,"z:\\projects\\boreal_swe_depletion\\data\\rep_subID.csv",sep=",",row.names=FALSE)

datRep <- dat.swe7[whichrep$rows,]



#pull out replicated data
repOut1 <- data.frame()
repOut2 <- data.frame()
repOut3 <- data.frame()
repmcmc <- mcmc.list()
repSum <- list()

for(i in 1:dim(chainDF)[1]){

	#read in output
	repOut1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\swerep_out.csv"))
	repOut2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\swerep_out.csv"))
	repOut3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\swerep_out.csv"))
	colnames(repOut1) <- paste0("rep_",seq(1,dim(repOut1)[2]))
	colnames(repOut2) <- paste0("rep_",seq(1,dim(repOut2)[2]))
	colnames(repOut3) <- paste0("rep_",seq(1,dim(repOut3)[2]))
	#turn into mcmc
	repOut1 <- as.mcmc(repOut1)
	repOut2 <- as.mcmc(repOut2)
	repOut3 <- as.mcmc(repOut3)
	repmcmc <- mcmc.list(repOut1,repOut2,repOut3)
	repSum[[i]] <- summary(repmcmc)
	}

#need to subset datrep
datRepL <- list()
for(i in 1:dim(chainDF)[1]){
		datRepL[[i]] <- datRep[datRep$gcID==chainDF$glc[i]&datRep$year==chainDF$year[i],]
}

#add rep data and ID	
repL <- list()
for(i in 1:dim(chainDF)[1]){

		repL[[i]] <- data.frame(repMean=repSum[[i]]$statistics[,1],
						datRepL[[i]])
						

}
repDF <- ldply(repL,data.frame)
colnames(midConv)[3] <- "gcID"
colnames(b0Conv)[3] <- "gcID"
colnames(b0Conv)[1] <- "b0conv"
colnames(b0Conv)[2] <- "pixID"
colnames(midConv)[2] <- "pixID"
#merge convergence info
repDF1 <- join(repDF,midConv,by=c("pixID","gcID","year"),type="left")
repDF2 <- join(repDF1,b0Conv,by=c("pixID","gcID","year"),type="left")



#only look at converged rep until those issues are solved
repDF3 <- repDF2[repDF2$conv==0&repDF2$b0conv==0,]

#double check this name
fit <- lm(repDF3$repMean~repDF3$sweN)
summary(fit)
plot(repDF3$sweN,repDF3$repMean, xlab="observed proportion of swe max",ylab="predicted proportion of swe max")
	
#check if there is something about the data in pixels not covnerging	


probMid <- midConv[midConv$conv==1,]
probB0 <- b0Conv[b0Conv$b0conv==1,]


dfConvM <- unique(data.frame(gcID=probMid$gcID,pixID=probMid$pixID)) 

dfProb <- join(probMid,pixJ3,by=c("gcID","pixID","year"),type="left")


plot(dat.swe5$jday[dat.swe5$pixID==probMid$pixID[1]&dat.swe5$year==probMid$year[1]&dat.swe5$gcID==probMid$gcID[1]],
	dat.swe5$sweN[dat.swe5$pixID==probMid$pixID[1]&dat.swe5$year==probMid$year[1]&dat.swe5$gcID==probMid$gcID[1]])
	
i=28
plot(dat.swe7$jday[dat.swe7$pixID==probMid$pixID[i]&dat.swe7$year==probMid$year[i]&dat.swe7$gcID==probMid$gcID[i]],
	dat.swe7$sweN[dat.swe7$pixID==probMid$pixID[i]&dat.swe7$year==probMid$year[i]&dat.swe7$gcID==probMid$gcID[i]])	
	
	
plot(dat.swe5$jday[dat.swe5$pixID==2&dat.swe5$year==2005&dat.swe5$gcID==3],
	dat.swe5$sweN[dat.swe5$pixID==2&dat.swe5$year==2005&dat.swe5$gcID==3])	
for(i in 1:dim(probMid)[1]){

points(dat.swe7$jday[dat.swe7$pixID==probMid$pixID[i]&dat.swe7$year==probMid$year[i]&dat.swe7$gcID==probMid$gcID[i]],
	dat.swe7$sweN[dat.swe7$pixID==probMid$pixID[i]&dat.swe7$year==probMid$year[i]&dat.swe7$gcID==probMid$gcID[i]], pch=19,type="l")	
}	
plot(dat.swe5$jday[dat.swe5$pixID==1&dat.swe5$year==2005&dat.swe5$gcID==1],
	dat.swe5$sweN[dat.swe5$pixID==1&dat.swe5$year==2005&dat.swe5$gcID==1])	
		
#make a map of cells

library(raster)
library(ncdf4)
library(rgdal)
library(gdalUtils)
#library(lubridate)

setwd("L:/data_repo/gis_data/")

# define the projection - EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 

########## SWE data from Mudryk ##########
# swe data in m from Mudryk
# daily mean SWE (GS2,MERRA2,Brown,Crocus)
# list the monthly files
swe.files <- list.files(pattern =".nc",path ="z:\\data_repo\\gis_data/swe_mudryk_blended\\",full.names=T)

# read one file in to use for reprojecting
pr <- raster(swe.files[1])
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')
plot(pr)	

#create unique ids with convergence info
dfConvM$convI <- rep(1,dim(dfConvM)[1])

cellInfo <- unique(data.frame(pixID=dat.swe5$pixID,gcID=dat.swe5$gcID,x.coord=dat.swe5$x.coord,y.coord=dat.swe5$y.coord))	

cellInfo <- join(cellInfo,dfConvM,by=c("pixID","gcID"),type="left")

cellIssue <- cellInfo[cellInfo$convI==1,]
points(cellInfo$x.coord,cellInfo$y.coord,pch=19,col="cornflowerblue", cex=.4)
points(cellIssue$x.coord,cellIssue$y.coord,pch=19,col="black",cex=.4)
legend("topleft",c("converged","convergence issue"),col=c("cornflowerblue","black"),pch=19)

#check data

plot(dat.swe5$jday[dat.swe5$pixID==1&dat.swe5$year==2000&dat.swe5$gcID==1],
	dat.swe5$sweN[dat.swe5$pixID==1&dat.swe5$year==2000&dat.swe5$gcID==1],
	xlab="Day of year",ylab="Normalized swe (proportion of maximum swe)")