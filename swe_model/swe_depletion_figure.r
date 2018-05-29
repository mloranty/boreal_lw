#######################################################
# Script for making figures                           #
# a model on a blended swe                            #
# product to look at differences in timing of         #
# swe depletion                                       #
#######################################################

library(plyr)
library(coda)
library(mcmcplots)

############################################
###  model input directory               ###
############################################
modDI <- "z:\\Projects\\boreal_swe_depletion\\model\\run3"
plotDI <- "z:\\Projects\\boreal_swe_depletion\\figures\\model\\run3"

############################################
###  read in swe data  and organize      ###
############################################
dat.swe <- read.csv("z:\\Projects\\boreal_swe_depletion\\data\\swe_depletion_model_data.csv")
dat.glc <- read.csv("z:\\Projects\\boreal_swe_depletion\\data\\glc50_table.csv")

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
dat.swe2 <- dat.swe1[dat.swe1$zone!=2&dat.swe1$zone!=3&dat.swe1$zone!=15&dat.swe1$zone!=16&dat.swe1$zone!=10&dat.swe1$zone!=14&
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
sweMaxF <- sweMax[sweMax$sweMax >=0.04,]
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

#join glc ID into dataframe
dat.swe4 <- join(dat.swe3,IDSglc, by="zone", type="left")
dat.swe4$t.airC <- dat.swe4$t.air-273.15

#get average temperature over a pixel, year

temp.py <- aggregate(dat.swe4$t.airC, by=list(dat.swe4$cell,dat.swe4$year,dat.swe4$zone), FUN="mean")
colnames(temp.py) <- c("cell","year","zone","temp")

#calculate average to center for each zone
temp.z <- aggregate(temp.py$temp,by=list(temp.py$zone),FUN="mean")
colnames(temp.z) <- c("zone","temp.zoneM")
#join zone mean back into temp.py
temp.py <- join(temp.py, temp.z, by="zone", type="left")
temp.py$tempCent <- temp.py$temp-temp.py$temp.zoneM

#join back into swe data
dat.swe5 <- join(dat.swe4,temp.py, by=c("cell","year","zone"), type="left")
############################################
###  read in model input                 ###
############################################
#list all output files

modFiles <- list.files(modDI)

#pull out parameter name
modParam <- character(0)
modGLC <- numeric(0)
modChainT <- character(0)
modChain <- numeric(0)
for(i in 1:length(modFiles)){
	modParam[i] <- strsplit(modFiles, "\\_")[[i]][1]
	modGLC[i] <- as.numeric(gsub("\\D","",strsplit(modFiles, "\\gc")[[i]][2]))
	modChainT[i] <- strsplit(modFiles, "\\gc")[[i]][1]
	modChain[i] <- as.numeric(gsub("\\D","",strsplit(modChainT,"\\chain")[[i]][2]))
	
}	

modelOut <- data.frame(parms=modParam,glc=modGLC,chain=modChain)
#subset to only look at glc 1 for now

#pull out each model chain
chain1 <- which(modelOut$chain==1&modelOut$glc==1)
chain2 <- which(modelOut$chain==2&modelOut$glc==1)
chain3 <- which(modelOut$chain==3&modelOut$glc==1)



#read in coda
chain1Out <- list()
chain2Out <- list()
chain3Out <- list()
for(i in 1:length(chain1)){
	chain1Out[[i]] <- read.csv(paste0(modDI,"\\",modFiles[chain1[i]]))
	colnames(chain1Out[[i]]) <- paste0(modelOut$parms[chain1[i]],modelOut$glc[chain1[i]])
	chain2Out[[i]] <- read.csv(paste0(modDI,"\\",modFiles[chain2[i]]))
	colnames(chain2Out[[i]]) <- paste0(modelOut$parms[chain2[i]],modelOut$glc[chain2[i]])
	chain3Out[[i]] <- read.csv(paste0(modDI,"\\",modFiles[chain3[i]]))
	colnames(chain3Out[[i]]) <- paste0(modelOut$parms[chain3[i]],modelOut$glc[chain3[i]])
}
############################################
###  turn into a coda object             ###
############################################
#organize coda
chain1all <- chain1Out[[1]]
chain2all <- chain2Out[[1]]
chain3all <- chain3Out[[1]]
for(i in 2:length(chain1)){
	chain1all <- cbind(chain1all,chain1Out[[i]])
	chain2all <- cbind(chain2all,chain2Out[[i]])
	chain3all <- cbind(chain3all,chain3Out[[i]])
	
}
chain1all <- as.mcmc(chain1all)
chain2all <- as.mcmc(chain2all)
chain3all <- as.mcmc(chain3all)

codaAll <- mcmc.list(chain1all,chain2all,chain3all)

modSum <- summary(codaAll)
datC <- data.frame(cbind(modSum$statistics,modSum$quantiles))
datC$glcID <- rep(1,times=9)
datC$parms <- gsub("\\d","",rownames(datC))

############################################
###  check history                       ###
############################################

mcmcplot(codaAll,dir=paste0(plotDI,"\\history") )


############################################
###  plot params                         ###
############################################
#make a simplier glc plotting name
IDSglc$name <- c("needleleaf deciduous","deciduous shrub","evergreen shrub","herbaceous",
					"needleleaf evergreen","mixed tree")
					
#normalized day
dat.swe4$dayN <- (dat.swe4$jday-32)/(182-32)
#swe curve function
sweC <- function(M, b,day,mid,base){
	M/(1+exp(b*(day-mid)))+base
}
#plotting info					
wd <- 50
hd <- 50

xl <- 0
xh <- 1
yl <- 0
yh <- 0.5
xplot <- seq(0,1,by=.01)

xL <- seq(32,182, by=10)
xS <- seq(0,1, length.out=length(xL))
yS <- seq(0,0.5,by=.05)
#dim(IDSglc)[1]
for(i in 1:1){
	jpeg(paste0(plotDI,"\\curves_",IDSglc$name[i],".jpeg"), width=2000,height=2000,quality=100)
	layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl,xh), ylim=c(yl,yh), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")
				
		points(dat.swe5$dayN,dat.swe5$swe,pch=19, col="grey75")

		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i],
							datC$Mean[rownames(datC)=="b01"][i],
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i],
							datC$Mean[rownames(datC)=="base01"][i]),type="l", lwd=6,
							col="goldenrod3")
							
							
		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*-1),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*-1),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*-1),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*-1)),type="l", lwd=6,
							col="royalblue")	

		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*-12),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*-12),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*-12),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*-12)),type="l", lwd=6,
							col="royalblue4")	
		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*1),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*1),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*1),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*1)),type="l", lwd=6,
							col="tomato3")	

			points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*14),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*14),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*14),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*14)),type="l", lwd=6,
							col="tomato4")							
		axis(1, xS,rep("", length(xS)), lwd.ticks=5)
		mtext(xL,at=xS, side=1, line=5, cex=2)
		axis(2,yS, rep("", length(yS)), lwd.ticks=5)
		mtext(yS,at=yS,side=2, line=5, cex=2, las=2)
	dev.off()
}					
					

hist(temp.py$tempCent[temp.py$zone==5])

############################################
###  plot data                           ###
############################################

##plot maximum swe value at beginning of season before melt
#get the maximum swe value by averaging well before melt

dat.sweMax <- dat.swe5[dat.swe5$jday<90,]
#aggregate by cell, year
sweMax <- aggregate(dat.sweMax$swe, by=list(dat.sweMax$cell,dat.sweMax$year,dat.sweMax$zone), FUN="mean")
colnames(sweMax) <- c("cell","year","zone","sweMax")

#join with temperature data
sweMax2 <- join(sweMax, temp.py, by=c("cell","year","zone"), type="left")
#join glc data
sweMax3 <- join(sweMax2, IDSglc, by=c("zone"),type="left")

#plotting info					
wd <- 50
hd <- 50

xl <- -20
xh <- 20
yl <- 0
yh <- 0.5
yearU <- list()
yearN <- numeric(0)
for(i in 1:dim(IDSglc)[1]){
	yearU[[i]] <- unique(sweMax3$year[sweMax3$gcID==i])
	yearN[i] <- length(yearU[[i]])
	jpeg(paste0(plotDI,"\\maxSWE\\maxSWE_",IDSglc$name[i],".jpeg"), width=2000,height=2000,quality=100)
	layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
		par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl,xh), ylim=c(yl,yh), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")
		for(j in 1:yearN[i]){		
			points(sweMax3$tempCent[sweMax3$gcID==i&sweMax3$year==yearU[[i]][j]],sweMax3$sweMax[sweMax3$gcID==i&sweMax3$year==yearU[[i]][j]],pch=19, col=j, cex=3)
		}
	dev.off()
}	
	