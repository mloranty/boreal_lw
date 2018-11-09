#######################################################
# Script for reading in a model run on a blended swe  #
# product to look at differences in timing of         #
# swe depletion                                       #
#######################################################

#create a vector with all of the parent directories
dirP <- c("z:\\projects\\boreal_swe_depletion\\model\\run9\\run1\\run1",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run2\\run2",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run3\\run3",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run4\\run4",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run5\\run5",
			"z:\\projects\\boreal_swe_depletion\\model\\run9\\run6\\run6")


#get a list of all of the files

dirA <- list()

for(i in 1:length(dirP)){
	dirA[[i]] <- list.dirs(paste0(dirP[i]), full.names=FALSE)

}