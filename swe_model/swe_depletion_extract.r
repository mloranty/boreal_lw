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
dirP <- c("z:\\projects\\boreal_swe_depletion\\model\\run10\\run1")

outD <- "z:\\projects\\boreal_swe_depletion\\model\\run10\\eval"
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
b0Sum <- list()

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
	b0Sum[[i]] <-cbind(summary(b0mcmc)$statistics,summary(b0mcmc)$quantiles)
	
}
	

midOut1 <- data.frame()
midOut2 <- data.frame()
midOut3 <- data.frame()
midmcmc <- mcmc.list()
midSum <- list()
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
	midSum[[i]] <-cbind(summary(midmcmc)$statistics,summary(midmcmc)$quantiles)
}