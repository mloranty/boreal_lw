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
rn <- 1


#linux =1 or windows =2
runOS <- 1
#linux data directory first option, windows second optioon
DDdir <- c("/mnt/g/projects/boreal_swe_depletion/data",
				"z:\\projects\\boreal_swe_depletion\\data")

#model directory
modDir <- "/home/hkropp/github/boreal_lw/swe_model/swe_depletion_model_code.stan"				

#output directory
outdir <- "/mnt/g/projects/boreal_swe_depletion/model/run2"


#######################################################
# read in and filter data                             #
#######################################################
#read in data files
if(runOS==1){
	dat.swe <- read.csv(paste0(DDdir[1], "/swe_depletion_model_data.csv"))
	dat.glc <- read.csv(paste0(DDdir[1], "/glc50_table.csv"))
}else{

	dat.swe <- read.csv(paste0(DDdir[2],"\\swe_depletion_model_data.csv"))
	dat.glc <- read.csv(paste0(DDdir[2], "/glc50_table.csv"))
}
print("finish reading in data")

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

print("finish data organize")
#######################################################
# set up model run                                    #
#######################################################



#inits for all possible runs					

inits2<-list(list(M=c(.1),base=c(.05),b=c(21),
	sig_swe=c(.05), mid=c(.4)))
inits3<-list(list(M=c(.08),base=c(.040),b=c(25)
		,sig_swe=c(.06), mid=c(.5)))		
inits1<-list(list(M=c(.07),base=c(.045),b=c(30),
			sig_swe=c(.1), mid=c(.6)))
			
			
for(i in 1:dim(IDSglc)[1]){		
			
	print(paste("start model run", i))			
	if(rn==1){		
	stan_model1 = stan(paste0(modDir), 
					data = list(Nobs=dim(dat.swe4[dat.swe4$gcID==i,])[1], swe=dat.swe4$swe[dat.swe4$gcID==i], 
				day=(dat.swe4$jday[dat.swe4$gcID==i]-32)/(182-32)),init=inits1,
				,chains=1, iter=3000)	
	print(paste("end model run",i))	
	out1<- extract(stan_model1)
	print(paste("extract variables",i))	
	write.table(out1$M, paste0(outdir,"/M_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$base, paste0(outdir,"/base_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$b, paste0(outdir,"/b_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$sig_swe, paste0(outdir,"/sig_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$mid, paste0(outdir,"/mid_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	print("end output")				
				
	}
	if(rn==2){			
	stan_model2 = stan(paste0(modDir), 
					data =  list(Nobs=dim(dat.swe4[dat.swe4$gcID==i,])[1], swe=dat.swe4$swe[dat.swe4$gcID==i], 
				day=(dat.swe4$jday[dat.swe4$gcID==i]-32)/(182-32)),init=inits2,
				,chains=1, iter=3000)	
	print(paste("end model run",i))
	out2<- extract(stan_model2)	
	print(paste("extract variables",i))	
	write.table(out2$M, paste0(outdir,"/M_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$base, paste0(outdir,"/base_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$b, paste0(outdir,"/b_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$sig_swe, paste0(outdir,"/sig_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$mid, paste0(outdir,"/mid_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	print(paste("end output",i))		

	}

	if(rn==3){	
	stan_model3 = stan(paste0(modDir), 
					data =  list(Nobs=dim(dat.swe4[sweDF3$gcID==i,])[1], swe=dat.swe4$swe[dat.swe4$gcID==i], 
				day=(dat.swe4$jday[dat.swe4$gcID==i]-32)/(182-32)),init=inits3,
				,chains=1, iter=3000)	
	print(paste("end model run",i))
	out3<- extract(stan_model3)	
	print(paste("extract variables",i))	
	write.table(out3$M, paste0(outdir,"/M_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$base, paste0(outdir,"/base_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$b, paste0(outdir,"/b_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$sig_swe, paste0(outdir,"/sig_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$mid, paste0(outdir,"/mid_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	print(paste("end output",i))		
				
	}
}

##this script ran on 806840 observations in the 6 classes in less than 2 hours