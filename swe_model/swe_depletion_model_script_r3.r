#######################################################
# Script for running a model on a blended swe         #
# product to look at differences in timing of         #
# swe depletion
#######################################################

#load libraries
library(coda)
library(plyr)
library(rstan)

#######################################################
# set up run info                                     #
#######################################################
#run number
rn <- 3


#linux =1 or windows =2
runOS <- 1
#linux data directory first option, windows second optioon
DDdir <- c("/mnt/g/projects/boreal_swe_depletion/data",
				"z:\\projects\\boreal_swe_depletion\\data")

#model directory
modDir <- "/home/hkropp/github/boreal_lw/swe_model/swe_depletion_model_code.stan"				

#output directory
outdir <- "/mnt/g/projects/boreal_swe_depletion/model/run5"


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
dat.swe <- dat.swe[dat.swe$year==2009,]




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

#normalize swe
#calculate percent
dat.swe5$sweP <- dat.swe5$swe/dat.swe5$sweMax

#round swe for 20% of peak to 1 and 20% of low to zero
dat.swe5$sweN <- ifelse(dat.swe5$sweP>=0.8,1,
					ifelse(dat.swe5$sweP<=0.2,0,dat.swe5$sweP))



#get unique pixel id in each glc for parameter id

pixID <- unique(data.frame(cell=dat.swe5$cell, gcID=dat.swe5$gcID))

#make sure pix doesn't change class
length(unique(dat.swe5$cell))
 dim(pixID)

#subset into each glc
pixList <- list()
pixGLC <- numeric(0)
for(i in 1:dim(IDSglc)[1]){
	pixList[[i]] <- pixID[pixID$gcID==i,]
	pixList[[i]]$pixID <- seq(1,dim(pixList[[i]])[1])
	pixGLC[i] <- dim(pixList[[i]])[1]
}

pixJ <- ldply(pixList,data.frame)

#join back into swe
dat.swe6 <- join(dat.swe5,pixJ, by=c("cell","gcID"), type="left")


print("finish data organize")
#######################################################
# set up model run                                    #
#######################################################




			
for(i in 1:dim(IDSglc)[1]){		
		
#inits for all possible runs					

inits2<-list(list(base0=rep(.05,pixGLC[i]),b0=c(21,pixGLC[i]),
	sig_swe=c(.05)))
inits3<-list(list(base0=c(.040,pixGLC[i]),b0=c(25,pixGLC[i])
		,sig_swe=c(.06)))		
inits1<-list(list(base0=c(.045,pixGLC[i]),b0=c(30,pixGLC[i]),
			sig_swe=c(.1)))
			

		
	print(paste("start model run", i))			
	if(rn==1){		
	stan_model1 = stan(paste0(modDir), 
					data = list(Nobs=dim(dat.swe6[dat.swe6$gcID==i,])[1], swe=dat.swe6$sweN[dat.swe6$gcID==i], 
				day=(dat.swe6$jday[dat.swe6$gcID==i]-32)/(182-32),pixID=dat.swe6$pixID,Npixel=pixGLC[i]),init=inits1,
				,chains=1, iter=3000)	
	print(paste("end model run",i))	
	out1<- extract(stan_model1)
	print(paste("extract variables",i))	
	write.table(out1$b0, paste0(outdir,"/b0_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$mid0, paste0(outdir,"/mid0_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	write.table(out1$sig_swe, paste0(outdir,"/sig_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	
	print("end output")				
				
	}
	if(rn==2){			
	stan_model2 = stan(paste0(modDir), 
					data =  list(Nobs=dim(dat.swe6[dat.swe6$gcID==i,])[1], swe=dat.swe6$sweN[dat.swe6$gcID==i], 
				day=(dat.swe6$jday[dat.swe6$gcID==i]-32)/(182-32),pixID=dat.swe6$pixID,Npixel=pixGLC[i]),chains=1, iter=3000)	
	print(paste("end model run",i))
	out2<- extract(stan_model2)	
	print(paste("extract variables",i))	

	write.table(out2$b0, paste0(outdir,"/b0_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$mid0, paste0(outdir,"/mid0_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	write.table(out2$sig_swe, paste0(outdir,"/sig_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	
	print(paste("end output",i))		

	}

	if(rn==3){	
	stan_model3 = stan(paste0(modDir), 
					data =  list(Nobs=dim(dat.swe6[dat.swe6$gcID==i,])[1], swe=dat.swe6$sweN[dat.swe6$gcID==i], 
				day=(dat.swe6$jday[dat.swe6$gcID==i]-32)/(182-32),pixID=dat.swe6$pixID,Npixel=pixGLC[i]),chains=1, iter=3000)	
	print(paste("end model run",i))
	out3<- extract(stan_model3)	
	print(paste("extract variables",i))	

	write.table(out3$b0, paste0(outdir,"/b0_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$mid0, paste0(outdir,"/mid0_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	write.table(out3$sig_swe, paste0(outdir,"/sig_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	print(paste("end output",i))		
				
	}
}

