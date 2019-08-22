##########################################################
#### Analysis of 9 years of swe depletion in boreal   ####
#### forests.                                         ####
#### outputs: daily swe data used in analysis: datSwe ####
#### and all swe before subsetting melt period: allSwe####
#### table of melt period rate and other info: cellSwe####
##########################################################

###############################################
### libraries                               ###
###############################################

library(plyr)

###############################################
### set up file paths                       ###
###############################################
swepath <- "z:\\data_repo\\gis_data"

#######################################################
# data info                                           #
#######################################################

#linux =1 or windows =2
runOS <- 2
#linux data directory first option, windows second optioon
DDdir <- c("/home/hkropp/boreal/data",
				"z:\\projects\\boreal_swe_depletion\\data")

#######################################################
# read in and filter data                             #
#######################################################
#read in data files
if(runOS==1){
	dat.swe <- read.csv(paste0(DDdir[1], "/swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[1], "/glc50_table.csv"))
	whichrep <- read.csv(paste0(DDdir[1],"/rep_subID_new.csv"))
	datExc <- read.csv(paste0(DDdir[1],"/prob_pix.csv"))

}else{

	dat.swe <- read.csv(paste0(DDdir[2],"\\swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[2], "/glc50_table.csv"))
	whichrep <- read.csv(paste0(DDdir[2],"\\rep_subID_new.csv"))
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
				
				
				
#######################################################
# calculate maximum and min swe info                  #
#######################################################

sweMax <- aggregate(dat.swe2$swe,by=list(dat.swe2$cell,dat.swe2$year),FUN="max")
colnames(sweMax) <- c("cell","year","sweMax")

#check the number of observations in data
nObs <- aggregate(dat.swe2$swe,by=list(dat.swe2$cell,dat.swe2$year),FUN="length")
colnames(nObs) <- c("cell","year","nObs")

#combine sweStats
sweSumm <- join(sweMax,nObs, by=c("cell","year"),type="full")


##########################
##### Filter point 3 #####
##########################				
#make sure there is a swe max more than 0.04				
#filter out areas with no snow extent in spring

sweSumm  <- sweSumm[sweSumm$sweMax > 0.04,]

##########################
##### Filter point 4 #####
##########################		
#only take sites with enough data for observation
sweSumm <- sweSumm[sweSumm$nObs >= 30,]
#calculate thresholds for rate calculation
sweSumm$maxThres <- 0.8*sweSumm$sweMax

#join back to swe to filter
dat.swe3 <- join(dat.swe2, sweSumm, by=c("cell","year"), type="inner")

#######################################################
# subset swe to only use the swe                      #
# after swe max is reached                            #
#######################################################
#swe melt period is defined as between last day of within 80% of max and fist day of less than 20% of max
#set up an index for whether the data is in swe max period
dat.swe3$maxID <- ifelse(dat.swe3$swe >= dat.swe3$maxThres,1,0)
#subset so only within max period
datMax <- dat.swe3[dat.swe3$maxID==1,]
#get onset of melt finding the last day of max swe
dayOnset <- aggregate(datMax$jday, by=list(datMax$cell,datMax$year), FUN="max",na.rm=TRUE)
colnames(dayOnset) <- c("cell","year","meltStart")
#join back into swe3
dat.swe3 <- join(dat.swe3,dayOnset,by=c("cell","year"),type="left")

#get max swe from lastmaxday
dat.sweMaxDay <- dat.swe3[dat.swe3$jday==dat.swe3$meltStart,]
#dataframe to join for swe on last day before melt
#since day 1 of melt is already the first day of melting occuring
MaximumSwe <- data.frame(sweStart=dat.sweMaxDay$swe,cell=dat.sweMaxDay$cell,year=dat.sweMaxDay$year)
 
dat.melt1 <- dat.swe3[dat.swe3$jday>=dat.swe3$meltStart,] 

sweMin <- aggregate(dat.melt1$swe ,by=list(dat.melt1$cell,dat.melt1$year),FUN="min")
colnames(sweMin) <- c("cell","year","sweMin")

#join back into sweSumm
sweSumm <- join(sweSumm,sweMin, by=c("cell","year"), type="full")
sweSumm <- join(sweSumm,MaximumSwe,by=c("cell","year"), type="full")

#calculate minimum threshold
sweSumm$minThres <- ifelse(0.2*sweSumm$sweMax < 0.01 | 0.2*sweSumm$sweMax < sweSumm$sweMin, sweSumm$sweMin, 0.2*sweSumm$sweMax)


#join back into swe melt
#mini df to join
minJoin <- data.frame(cell=sweSumm$cell,year=sweSumm$year,minThres=sweSumm$minThres)

dat.melt1 <- join(dat.melt1,minJoin,by=c("cell","year"),type="left")

#set up threshold for whether the data is in the swe min period
dat.melt1$minID <- ifelse(dat.melt1$swe <= dat.melt1$minThres,1,0) 

#get last day of melt finding the first day of swe min
datMin <- dat.melt1[dat.melt1$minID==1,]
#get end of melt finding the first day of min swe
dayEnd <- aggregate(datMin$jday, by=list(datMin$cell,datMin$year), FUN="min",na.rm=TRUE)
colnames(dayEnd) <- c("cell","year","meltEnd")


#join back into datmelt
dat.melt2 <- join(dat.melt1,dayEnd, by=c("cell","year"),type="left")

#get swe on last day of melt
dat.sweMinDay <- dat.melt2[dat.melt2$jday==dat.melt2$meltEnd,]
#join actual swe in
sweMinDayJ <- data.frame(sweEnd=dat.sweMinDay$swe,cell=dat.sweMinDay$cell,year=dat.sweMinDay$year)

#join back into sweSumm
sweSumm <- join(sweSumm,dayOnset, by=c("cell","year"), type="left")
sweSumm <- join(sweSumm,dayEnd, by=c("cell","year"), type="left")
sweSumm <- join(sweSumm,sweMinDayJ, by=c("cell","year"), type="left")
#calculate the length of the melt
sweSumm$meltDays <- sweSumm$meltEnd-sweSumm$meltStart

#check melt days
if(nrow(sweSumm[is.na(sweSumm$meltDays),])==0){
	print("no missing melt days due to error")
	}else{
	print("Error missing melt days due to error")
	}
#check negative days
if(nrow(sweSumm[sweSumm$meltDays<0,])==0){
	print("no negative days in melt period")
}else{
	print("Error: negative days in melt period")
	}

##########################
##### Filter point 5 #####
##########################	
# a melt period of less than 5 days is  not appropriate for this analysis
sweSumm <- sweSumm[sweSumm$meltDays > 5,]

#calculate rate of swe depletion during melt period
sweSumm$meltRate <- (sweSumm$sweStart-sweSumm$sweEnd)/(	sweSumm$meltStart-sweSumm$meltEnd)
sweSumm$meltRateCM <- ((sweSumm$sweStart-sweSumm$sweEnd)*100)/(	sweSumm$meltStart-sweSumm$meltEnd)

#check no positive melt rates
if(nrow(sweSumm[sweSumm$meltRate>0,])==0){
	print("no positive melt rates")
	}else{
		print("Error: positive melt rates")
	}


#finish sub-setting melt
dat.melt3 <- dat.melt2[dat.melt2$jday<=dat.melt2$meltEnd,]

#organize glc IDS
IDSglc <- unique(data.frame(zone=dat.melt3$zone))
IDSglc <- join(IDSglc, dat.glc[,1:2], by="zone", type="left")
IDSglc$gcID <- seq(1,dim(IDSglc)[1])

#add glc, tree cover into sweSumm
cellInfo <- unique(data.frame(cell=dat.swe3$cell,year=dat.swe3$year, vcf=dat.swe3$vcf, zone=dat.swe3$zone))
#join into summ table. This will also be final table
cellSwe <- join(sweSumm,cellInfo, by=c("cell","year"),type="left")

#calculate average surface temperature during the melt period
surfTemp <- aggregate(dat.melt3$t.air-273.15, by=list(dat.melt3$cell,dat.melt3$year), FUN="mean",na.rm=TRUE)
colnames(surfTemp) <- c("cell","year","tair")

cellSwe <- join(cellSwe,surfTemp,by=c("cell","year"),type="left")

#check no missing tair
if(nrow(cellSwe[is.na(cellSwe$tair),])==0){
	print("no missing surface temp due to error")
	}else{
	print("Error: missing surface temp")
	}

#organize data names
#swe during melt period
datSwe <- dat.melt3
#all swe
sweAll <- dat.swe3 


print("finish data organize")

rm(list=setdiff(ls(), c("datSwe","IDSglc","sweAll", "cellSwe")))