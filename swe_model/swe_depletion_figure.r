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
modDI <- "z:\\Projects\\boreal_swe_depletion\\model\\run1\\run1"

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

mcmcplot(codaAll,dir=paste0(modDI,"\\history") )


############################################
###  plot params                         ###
############################################

