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

#create a vector with all of the parent directories
dirP <- c("z:\\projects\\boreal_swe_depletion\\model\\run9\\run1\\run1",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run2\\run2",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run3\\run3",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run4\\run4",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run5\\run5",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run6\\run6")

outD <- "z:\\projects\\boreal_swe_depletion\\model\\run9\\eval"
#get a list of all of the files

dirA <- list()
dimA <- numeric(0)
datA <- list()
for(i in 1:length(dirP)){
	dirA[[i]] <- list.dirs(paste0(dirP[i]), full.names=FALSE)
	#remove parent directory
	dirA[[i]] <- dirA[[i]][-1]
	dimA[i] <- length(dirA[[i]])
	datA[[i]] <- data.frame(files=dirA[[i]],dirP=rep(dirP[i],dimA[i]))
	
	
}

datA <- ldply(datA,data.frame)
#create a 

#pull out identifying info

for(i in 1:dim(datA)[1]){
		#vegetation land class
		datA$glc[i] <- as.numeric(strsplit(as.character(datA$files[i]),"_")[[1]][2])
		#year
		datA$year[i] <- as.numeric(strsplit(as.character(datA$files[i]),"_")[[1]][4])
		#chain
		datA$chain[i] <- as.numeric(strsplit(as.character(datA$files[i]),"_")[[1]][6])
}	

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
	#dir.create(paste0(outD,"\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
	#mcmcplot(b0mcmc,dir=paste0(outD,"\\glc",chainDF$glc[i],"_year",chainDF$year[i]))
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


#######################################################
# read in and filter data                             #
#######################################################
#read in data files
if(runOS==1){
	dat.swe <- read.csv(paste0(DDdir[1], "/swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[1], "/glc50_table.csv"))

}else{

	dat.swe <- read.csv(paste0(DDdir[2],"\\swe_depletion_model_data_vcf_no_topo.csv"))
	dat.glc <- read.csv(paste0(DDdir[2], "/glc50_table.csv"))

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
#join mod run with gcID
colnames(chainDF)[3] <- "gcID" 
gcYearID <- join(gcYearID,chainDF, by=c("gcID","year"), type="left")

datalist <- list()
rep1 <- matrix(rep(NA,6000*1000),ncol=6000)
#for(i in 1:dim(gcYearID)[1]){
ptm <- proc.time()
proc.time() - ptm
for(i in 1:1){
		datalist[[i]] <- list(Nobs=dim(dat.swe5[dat.swe5$gcID==gcYearID$gcID[i]&dat.swe5$year==gcYearID$year[i],])[1],
				swe=dat.swe5$sweN[dat.swe5$gcID==gcYearID$gcID[i]&dat.swe5$year==gcYearID$year[i]], 
				day=(dat.swe5$jday[dat.swe5$gcID==gcYearID$gcID[i]&dat.swe5$year==gcYearID$year[i]]-32)/(182-32),
				pixID=dat.swe5$pixID[dat.swe5$gcID==gcYearID$gcID[i]&dat.swe5$year==gcYearID$year[i]],
				Npixel=dim(pixList[[i]])[1])
				
			#read in output
		b0Out1 <- read.csv(paste0(gcYearID$dirP1[i],"\\",gcYearID$files1[i],"\\b0_out.csv"))
		midOut1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\mid0_out.csv"))
		sigOut1 <- read.csv(paste0(chainDF$dirP1[i],"\\",chainDF$files1[i],"\\sig_out.csv"))
		b0Out2 <- read.csv(paste0(gcYearID$dirP2[i],"\\",gcYearID$files2[i],"\\b0_out.csv"))
		midOut2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\mid0_out.csv"))
		sigOut2 <- read.csv(paste0(chainDF$dirP2[i],"\\",chainDF$files2[i],"\\sig_out.csv"))
		b0Out3 <- read.csv(paste0(gcYearID$dirP3[i],"\\",gcYearID$files3[i],"\\b0_out.csv"))
		midOut3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\mid0_out.csv"))
		sigOut3 <- read.csv(paste0(chainDF$dirP3[i],"\\",chainDF$files3[i],"\\sig_out.csv"))
		b0Out <- rbind(b0Out1,b0Out2,b0Out3)
		midOut <- rbind(midOut1,midOut2,midOut3)
		sigOut <- rbind(sigOut1,sigOut2,sigOut3)
		#calculate rep	
		ptm <- proc.time()
		repl <- matrix(rep(NA,6000*1000),ncol=6000)
		for(j in 1:dim(b0Out)[1]){
			#for(m in 1:dim(datalist[[i]])[1]){
			for(m in 1:1000){
				rep1[m,j] <- rnorm(1,1/(1+exp(b0Out[j,datalist[[i]]$pixID[m]]*(datalist[[i]]$day[m]-midOut[j,datalist[[i]]$pixID[m]]))), sigOut[j,1])
			}
		}
}
proc.time() - ptm
rnorm(1,1/(1+exp(b0Out1[1,datalist[[1]]$pixID[1]]*(datalist[[1]]$day[1]-midOut1[1,datalist[[1]]$pixID[1]]))), 1)
#normal(1/(1+exp(b0[pixID[i]]*(day[i]-mid0[pixID[i]]))), sig_swe)