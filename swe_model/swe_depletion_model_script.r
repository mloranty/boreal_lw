#!/shared/R-3.4.3/bin/Rscript
#######################################################
# Script for running a model on a blended swe         #
# product to look at differences in timing of         #
# swe depletion
#######################################################

#load libraries
#library(coda, lib="/home/hkropp/Rtemp/R3.4.4")
#library(plyr, lib="/home/hkropp/Rtemp/R3.4.4")
library(coda)
library(plyr)
library(rstan)
library(snow)
library(snowfall)
#library(snow, lib="/home/hkropp/Rtemp/R3.4.4")
#library(snowfall, lib="/home/hkropp/Rtemp/R3.4.4")

#######################################################
# set up run info                                     #
#######################################################
#run number
rn <- 1
#chain number
chain <- 1

#linux =1 or windows =2
runOS <- 1
#linux data directory first option, windows second optioon
DDdir <- c("/home/hkropp/boreal/data",
				"z:\\projects\\boreal_swe_depletion\\data")

#model directory
modDir <- "/home/hkropp/github/boreal_lw/swe_model/swe_depletion_model_code.stan"				

#output directory
outdir <- "/home/hkropp/boreal/model/run1/"


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
	whichrep <- read.csv(paste0(DDdir[1],"\\rep_subID.csv"))
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

#######################################################
# set up model run                                    #
#######################################################
#need to subset data for each glc and year
datalist <- list()
if(rn==1){
	for(i in 1:30){

		datalist[[i]] <- list(Nobs=dim(dat.swe7[dat.swe7$gcID==gcYearID$gcID[i]&dat.swe7$year==gcYearID$year[i],])[1],
				swe=dat.swe7$sweN[dat.swe7$gcID==gcYearID$gcID[i]&dat.swe7$year==gcYearID$year[i]], 
				day=(dat.swe7$jday[dat.swe7$gcID==gcYearID$gcID[i]&dat.swe7$year==gcYearID$year[i]]-32)/(182-32),
				pixID=dat.swe7$pixID[dat.swe7$gcID==gcYearID$gcID[i]&dat.swe7$year==gcYearID$year[i]],
				Npixel=dim(pixList[[i]])[1], 
				Rday=(datRep$jday[datRep$gcID==gcYearID$gcID[i]&datRep$year==gcYearID$year[i]]-32)/(182-32),
				Nrep=dim(datRep[datRep$gcID==gcYearID$gcID[i]&datRep$year==gcYearID$year[i],])[1],
				RpixID=datRep$pixID[datRep$gcID==gcYearID$gcID[i]&datRep$year==gcYearID$year[i]])
	}
}

if(rn==2){
	for(i in 1:20){
				datalist[[i]] <- list(Nobs=dim(dat.swe7[dat.swe7$gcID==gcYearID$gcID[i+30]&dat.swe7$year==gcYearID$year[i+30],])[1],
				swe=dat.swe7$sweN[dat.swe7$gcID==gcYearID$gcID[i+30]&dat.swe7$year==gcYearID$year[i+30]], 
				day=(dat.swe7$jday[dat.swe7$gcID==gcYearID$gcID[i+30]&dat.swe7$year==gcYearID$year[i+30]]-32)/(182-32),
				pixID=dat.swe7$pixID[dat.swe7$gcID==gcYearID$gcID[i+30]&dat.swe7$year==gcYearID$year[i+30]],
				Npixel=dim(pixList[[i+30]])[1],  
				Rday=(datRep$jday[datRep$gcID==gcYearID$gcID[i+30]&datRep$year==gcYearID$year[i+30]]-32)/(182-32),
				Nrep=dim(datRep[datRep$gcID==gcYearID$gcID[i+30]&datRep$year==gcYearID$year[i+30],])[1],
				RpixID=datRep$pixID[datRep$gcID==gcYearID$gcID[i+30]&datRep$year==gcYearID$year[i+30]])
	}
}



#set up inits 
if(chain==1){
	inits <- list(list(mu_b0=50,sig_b0=10,mu_mid=.5,sig_mid=.05))
}

if(chain==2){
	inits <- list(list(mu_b0=25,sig_b0=5,mu_mid=.7,sig_mid=.1))
}

if(chain==3){
	inits <- list(list(mu_b0=75,sig_b0=15,mu_mid=.6,sig_mid=.2))
}

# set the number of CPUs to be 32 for a node on the cluster

if(rn==1){
	sfInit(parallel=TRUE, cpus=30)

}


if(rn==2){
sfInit(parallel=TRUE, cpus=20)
}


#assign rstan to each CPU
sfLibrary(rstan)

#create unique filepath for each model run
dirList <- list()
if(rn==1){
	for(i in 1:30){
		dirList[[i]] <- paste0(outdir,"gcID_",gcYearID$gcID[i],"_year_",gcYearID$year[i],"_chain_",chain,"_run")
		dir.create(dirList[[i]])
	}


}

if(rn==2){
	for(i in 1:20){
	dirList[[i]] <- paste0(outdir,"gcID_",gcYearID$gcID[i+30],"_year_",gcYearID$year[i+30],"_chain_",chain,"_run")
	dir.create(dirList[[i]])
	}
}

#set up stan runction	
parallel.stan <- function(X,dataL,init,outDIR){
			stan_model1 = stan(file="/home/hkropp/github/boreal_lw/swe_model/swe_depletion_model_code.stan", 
					data = dataL[[X]], init=init,
				,chains=1, warmup=8000, iter=10000)
			out1= extract(stan_model1)
			write.table(out1$sig_swe, paste0(outDIR[[X]],"/sig_out.csv"), sep=",")
			write.table(out1$b0, paste0(outDIR[[X]],"/b0_out.csv"), sep=",")
			write.table(out1$mid0, paste0(outDIR[[X]],"/mid0_out.csv"), sep=",")
			write.table(out1$muB0, paste0(outDIR[[X]],"/muB0_out.csv"), sep=",")
			write.table(out1$sigB0, paste0(outDIR[[X]],"/sigB0_out.csv"), sep=",")
			write.table(out1$muMid, paste0(outDIR[[X]],"/muMid_out.csv"), sep=",")
			write.table(out1$sigMid, paste0(outDIR[[X]],"/sigMid_out.csv"), sep=",")
			write.table(out1$swerep, paste0(outDIR[[X]],"/swerep_out.csv"), sep=",")
}

if(rn==1){
	sfLapply(1:30, parallel.stan,dataL=datalist,init=inits,outDIR=dirList)	
	
}	

if(rn==2){
	sfLapply(1:20, parallel.stan,dataL=datalist,init=inits,outDIR=dirList)			
}	

sfStop()
