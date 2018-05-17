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
modDI <- "z:\\Projects\\boreal_swe_depletion\\model\\run2"
plotDI <- "z:\\Projects\\boreal_swe_depletion\\figures\\model\\run2"

############################################
###  read in swe data  and organize      ###
############################################
dat.swe <- read.csv("z:\\Projects\\boreal_swe_depletion\\data\\swe_depletion_model_data.csv")
dat.glc <- read.csv("z:\\Projects\\boreal_swe_depletion\\data\\swe_depletion_model_data.csv")

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

#pull out each model chain
chain1 <- which(modelOut$chain==1)
chain2 <- which(modelOut$chain==2)
chain3 <- which(modelOut$chain==3)
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
############################################
###  check history                       ###
############################################

mcmcplot(codaAll,dir=paste0(plotDI,"\\history") )


############################################
###  plot params                         ###
############################################


