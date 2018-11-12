#######################################################
# Script for reading in a model run on a blended swe  #
# product to look at differences in timing of         #
# swe depletion                                       #
#######################################################
library(plyr)



######################################
#### file name information        ####
######################################


#create a vector with all of the parent directories
dirP <- c("z:\\projects\\boreal_swe_depletion\\model\\run9\\run1\\run1",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run2\\run2",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run3\\run3",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run4\\run4",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run5\\run5",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run6\\run6")


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
