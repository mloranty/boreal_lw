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
dirP <- c("z:\\projects\\boreal_swe_depletion\\model\\run14\\run1")

outD <- "z:\\projects\\boreal_swe_depletion\\model\\run14\\eval"
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
whichrep <- read.csv("z:\\projects\\boreal_swe_depletion\\data\\rep_subID_new.csv")
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
	dir.create(paste0(outD,"\\mid\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
	mcmcplot(midmcmc,dir=paste0(outD,"\\mid\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
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
		datRepL[[i]] <- whichrep[whichrep$gcID==chainDF$glc[i]&whichrep$year==chainDF$year[i],]
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
repDF1 <- join(repDF,midConv,by=c("gcID","year"),type="left")
repDF2 <- join(repDF1,b0Conv,by=c("gcID","year"),type="left")



#only look at converged rep until those issues are solved
#repDF3 <- repDF2[repDF2$conv==0&repDF2$b0conv==0,]

#double check this name
fit <- lm(repDF2$repMean~repDF2$sweN)
summary(fit)
plot(repDF2$sweN,repDF2$repMean, xlab="observed proportion of swe max",ylab="predicted proportion of swe max")
abline(0,1, col="red",lwd=2)	
abline(fit, col="cornflowerblue",lwd=2,lty=3)
legend("topleft",c("1:1 line","regression line"),col=c("red","cornflowerblue"),lwd=2,lty=c(1,2), bty="n")
#check if there is something about the data in pixels not covnerging	


probMid <- midConv[midConv$conv==1,]
probB0 <- b0Conv[b0Conv$conv==1,]

probAll <- join(probMid,probB0, by=c("pixID","gcID","year"),type="full")

probOut <- data.frame(pixID=probAll$pixID, gcID=probAll$gcID,year=probAll$year)

#write.table (probOut, "z:\\projects\\boreal_swe_depletion\\data\\prob_pix.csv",sep=",", row.names=FALSE)
#note three sites show up as marginally not making the gelman metric but they have output chains that converged
#but have very long tails, likely making the variance large in the metric. Including these in analysis

#write ID to file
		
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