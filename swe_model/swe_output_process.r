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
dirP <- c("z:\\projects\\boreal_swe_depletion\\model\\run13\\run1")

outD <- "z:\\projects\\boreal_swe_depletion\\model\\run13\\eval"
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
b0Summ <- list()
b0Stat <- list()
for(i in 1:dim(chainDF)[1]){

	#read in output
	b0Out1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\b0_out.csv"))
	b0Out2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\b0_out.csv"))
	b0Out3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\b0_out.csv"))
	colnames(b0Out1) <- paste0("b0_",seq(1,dim(b0Out1)[2]))
	colnames(b0Out2) <- paste0("b0_",seq(1,dim(b0Out2)[2]))
	colnames(b0Out3) <- paste0("b0_",seq(1,dim(b0Out3)[2]))
	
	#turn into mcmc
	#and convert to days
	b0Out1 <- as.mcmc(apply(b0Out1,c(1,2),function(X){X/(182-32)}))
	b0Out2 <- as.mcmc(apply(b0Out2,c(1,2),function(X){X/(182-32)}))
	b0Out3 <- as.mcmc(apply(b0Out3,c(1,2),function(X){X/(182-32)}))
	b0mcmc <- mcmc.list(b0Out1,b0Out2,b0Out3)
	b0Summ[[i]] <- summary(b0mcmc)
	b0Stat[[i]] <- data.frame(b0Summ[[i]]$statistics,b0Summ[[i]]$quantiles,gcID=rep(chainDF$glc[i],dim(b0Summ[[i]]$statistics)[1]),
					year=rep(chainDF$year[i],dim(b0Summ[[i]]$statistics)[1]),pixID=seq(1,dim(b0Summ[[i]]$statistics)[1]))
}
	

midOut1 <- data.frame()
midOut2 <- data.frame()
midOut3 <- data.frame()
midmcmc <- mcmc.list()
mid.diag <- list()
midSumm <- list()
midStat <- list()
for(i in 1:dim(chainDF)[1]){

	#read in output
	midOut1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\mid0_out.csv"))
	midOut2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\mid0_out.csv"))
	midOut3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\mid0_out.csv"))
	colnames(midOut1) <- paste0("mid0_",seq(1,dim(midOut1)[2]))
	colnames(midOut2) <- paste0("mid0_",seq(1,dim(midOut2)[2]))
	colnames(midOut3) <- paste0("mid0_",seq(1,dim(midOut3)[2]))

	#turn into mcmc
	#convert to doy
	midOut1 <- as.mcmc(apply(midOut1,c(1,2),function(X){(X*(182-32))+32}))
	midOut2 <- as.mcmc(apply(midOut2,c(1,2),function(X){(X*(182-32))+32}))
	midOut3 <- as.mcmc(apply(midOut3,c(1,2),function(X){(X*(182-32))+32}))
	midmcmc <- mcmc.list(midOut1,midOut2,midOut3)
	midSumm[[i]] <- summary(midmcmc)
	midStat[[i]] <- data.frame(midSumm[[i]]$statistics,midSumm[[i]]$quantiles,gcID=rep(chainDF$glc[i],dim(midSumm[[i]]$statistics)[1]),
					year=rep(chainDF$year[i],dim(midSumm[[i]]$statistics)[1]),pixID=seq(1,dim(midSumm[[i]]$statistics)[1]))
}

midOut <- ldply(midStat, data.frame)
b0Out <- ldply(b0Stat, data.frame)

###
#read in means
mumidOut1 <- data.frame()
mumidOut2 <- data.frame()
mumidOut3 <- data.frame()
mumidmcmc <- mcmc.list()
mumidSumm <- list()
mumidStat <- list()
for(i in 1:dim(chainDF)[1]){
	#read in output
	mumidOut1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\muMid_out.csv"))
	mumidOut2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\muMid_out.csv"))
	mumidOut3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\muMid_out.csv"))
	colnames(mumidOut1) <- paste0("mumid")
	colnames(mumidOut2) <- paste0("mumid")
	colnames(mumidOut3) <- paste0("mumid")
	#turn into mcmc
	mumidOut1 <- as.mcmc(apply(mumidOut1,c(1,2),function(X){(X*(182-32))+32}))
	mumidOut2 <- as.mcmc(apply(mumidOut2,c(1,2),function(X){(X*(182-32))+32}))
	mumidOut3 <- as.mcmc(apply(mumidOut3,c(1,2),function(X){(X*(182-32))+32}))
	mumidmcmc <- mcmc.list(mumidOut1,mumidOut2,mumidOut3)
	mumidSumm[[i]] <- summary(mumidmcmc)
	mumidStat[[i]] <- data.frame(Mean=mumidSumm[[i]]$statistics[1],
									SD=mumidSumm[[i]]$statistics[2],
									p2.5=mumidSumm[[i]]$quantiles[1],
									p97.5=mumidSumm[[i]]$quantiles[5],gcID=chainDF$glc[i],
					year=chainDF$year[i])
}

muMidOut <- ldply(mumidStat,data.frame)

mub0Out1 <- data.frame()
mub0Out2 <- data.frame()
mub0Out3 <- data.frame()
mub0mcmc <- mcmc.list()
mub0.diag <- list()
mub0Summ <- list()
mub0Stat <- list()
for(i in 1:dim(chainDF)[1]){

	#read in output
	mub0Out1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\muB0_out.csv"))
	mub0Out2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\muB0_out.csv"))
	mub0Out3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\muB0_out.csv"))
	colnames(mub0Out1) <- paste0("mub0")
	colnames(mub0Out2) <- paste0("mub0")
	colnames(mub0Out3) <- paste0("mub0")
	#turn into mcmc
	mub0Out1 <- as.mcmc(apply(mub0Out1,c(1,2),function(X){X/(182-32)}))
	mub0Out2 <- as.mcmc(apply(mub0Out2,c(1,2),function(X){X/(182-32)}))
	mub0Out3 <- as.mcmc(apply(mub0Out3,c(1,2),function(X){X/(182-32)}))
	mub0mcmc <- mcmc.list(mub0Out1,mub0Out2,mub0Out3)
	mub0Summ[[i]] <- summary(mub0mcmc)
	mub0Stat[[i]] <- data.frame(Mean=mub0Summ[[i]]$statistics[1],
									SD=mub0Summ[[i]]$statistics[2],
									p2.5=mub0Summ[[i]]$quantiles[1],
									p97.5=mub0Summ[[i]]$quantiles[5],
							gcID=chainDF$glc[i],
					year=chainDF$year[i])
}
muB0Out <- ldply(mub0Stat,data.frame)
#######################################################
# read in and filter data                             #
#######################################################
#read in data files
if(runOS==1){
	dat.swe <- read.csv(paste0(DDdir[1], "/swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[1], "/glc50_table.csv"))
	whichrep <- read.csv(paste0(DDdir[1],"/rep_subID.csv"))
	datExc <- read.csv(paste0(DDdir[1],"/prob_pix.csv"))

}else{

	dat.swe <- read.csv(paste0(DDdir[2],"\\swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[2], "/glc50_table.csv"))
	whichrep <- read.csv(paste0(DDdir[2],"\\rep_subID.csv"))
	datExc <- read.csv(paste0(DDdir[2],"\\prob_pix.csv"))
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


################################################
####omit data after zero is reached         ####
################################################
minTemp <- numeric(0)
minN <- numeric(0)
#get the final swe max time
for(i in 1:dim(gcYearID)[1]){
	minTemp <- numeric(0)
	minN <- numeric(0)
	for(j in 1:dim(pixList[[i]])[1]){
		if(length(which(dat.swe7$pixID==pixList[[i]]$pixID[j]&dat.swe7$year==pixList[[i]]$year[j]& dat.swe7$gcID==pixList[[i]]$gcID[j]&dat.swe7$sweN==0))!=0){
			minTemp <- which(dat.swe7$pixID==pixList[[i]]$pixID[j]&dat.swe7$year==pixList[[i]]$year[j]& dat.swe7$gcID==pixList[[i]]$gcID[j]&dat.swe7$sweN==0) 
		
			minN[j] <- head(minTemp, n=1)
		}else{minN[j] <- NA}	
	}
	pixList[[i]]$finalMin <- minN
}

pixJ3 <- ldply(pixList,data.frame)

pixJ3$dayMin <- dat.swe7$jday[pixJ3$finalMin]
pixJ3$dayMax <- dat.swe5$jday[pixJ2$finalMax]

dat.swe8 <- join(dat.swe7,pixJ3, by=c("cell","year","gcID","pixID","gcYearID"), type="left")
dat.swe8$filterID <- ifelse(dat.swe8$jday<=dat.swe8$dayMin|is.na(dat.swe8$dayMin),1,0)

dat.swe9 <- dat.swe8[dat.swe8$filterID==1,]

################################################
####omit sites with a lot of missing data   ####
################################################
#get the maximum day of year
nDOY <- numeric(0)
#get the final swe max time
for(i in 1:dim(gcYearID)[1]){
	nDOY <- numeric(0)
	for(j in 1:dim(pixList[[i]])[1]){
		nDOY[j] <- length(dat.swe9$jday[dat.swe9$pixID==pixList[[i]]$pixID[j]&dat.swe9$year==pixList[[i]]$year[j]&dat.swe9$gcID==pixList[[i]]$gcID[j]])
		}
		pixList[[i]]$doyN <- nDOY
	}
pixJ4 <- ldply(pixList,data.frame)
pixJ4$dayMax <- dat.swe5$jday[pixJ2$finalMax]

dat.swe10 <- join(dat.swe9,pixJ4, by=c("cell","year","gcID","pixID","gcYearID"), type="left")


dat.swe11<- dat.swe10[dat.swe10$doyN>10,]	

datSwe <- dat.swe11
#clear out files that aren't needed
rm(list=setdiff(ls(), c("datSwe","b0Out","midOut","muB0Out","muMidOut","IDSglc")))
