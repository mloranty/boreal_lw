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
outdir <- "/mnt/g/projects/boreal_swe_depletion/model/run4"


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

#just run the model for a few years
dat.swe5 <- dat.swe5[dat.swe5$year>=2013,]


print("finish data organize")
#######################################################
# set up model run                                    #
#######################################################



#inits for all possible runs					

inits2<-list(list(M0=c(.1),base0=c(.05),b0=c(21),
	sig_swe=c(.05), mid0=c(.4), M1=c(0.001),base1=c(0.001),b1=c(0.0001),mid1=c(0.001),M2=c(0.001),base2=c(0.001),b2=c(0.0001),mid2=c(0.001)))
inits3<-list(list(M0=c(.08),base0=c(.040),b0=c(25)
		,sig_swe=c(.06), mid0=c(.5), M1=c(0.002),base1=c(0.002),b1=c(0.0002),mid1=c(0.002),
		M2=c(0.002),base2=c(0.002),b2=c(0.0002),mid2=c(0.002)))		
inits1<-list(list(M0=c(.07),base0=c(.045),b0=c(30),
			sig_swe=c(.1), mid0=c(.6), M1=c(0.003),base1=c(0.003),b1=c(0.0003),mid1=c(0.003),
			 M2=c(0.003),base2=c(0.003),b2=c(0.0003),mid2=c(0.003)))
			
			
for(i in 1:dim(IDSglc)[1]){		
			
	print(paste("start model run", i))			
	if(rn==1){		
	stan_model1 = stan(paste0(modDir), 
					data = list(Nobs=dim(dat.swe5[dat.swe5$gcID==i,])[1], swe=dat.swe5$swe[dat.swe5$gcID==i], 
				day=(dat.swe5$jday[dat.swe5$gcID==i]-32)/(182-32),tempC=dat.swe5$tempCent[dat.swe5$gcID==i],
				vcf=dat.swe5$vcf[dat.swe5$gcID==i]),init=inits1,
				,chains=1, iter=3000)	
	print(paste("end model run",i))	
	out1<- extract(stan_model1)
	print(paste("extract variables",i))	
	write.table(out1$M0, paste0(outdir,"/M0_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$base0, paste0(outdir,"/base0_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$b0, paste0(outdir,"/b0_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$mid0, paste0(outdir,"/mid0_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	write.table(out1$sig_swe, paste0(outdir,"/sig_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$M1, paste0(outdir,"/M1_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$base1, paste0(outdir,"/base1_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$b1, paste0(outdir,"/b1_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$mid1, paste0(outdir,"/mid1_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$M2, paste0(outdir,"/M2_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$base2, paste0(outdir,"/base2_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$b2, paste0(outdir,"/b2_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$mid2, paste0(outdir,"/mid2_out_chain1_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	
	print("end output")				
				
	}
	if(rn==2){			
	stan_model2 = stan(paste0(modDir), 
					data =  list(Nobs=dim(dat.swe5[dat.swe5$gcID==i,])[1], swe=dat.swe5$swe[dat.swe5$gcID==i], 
				day=(dat.swe5$jday[dat.swe5$gcID==i]-32)/(182-32),tempC=dat.swe5$tempCent[dat.swe5$gcID==i],vcf=dat.swe5$vcf[dat.swe5$gcID==i]),init=inits2,
				,chains=1, iter=3000)	
	print(paste("end model run",i))
	out2<- extract(stan_model2)	
	print(paste("extract variables",i))	
	write.table(out2$M0, paste0(outdir,"/M0_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$base0, paste0(outdir,"/base0_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$b0, paste0(outdir,"/b0_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$mid0, paste0(outdir,"/mid0_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	write.table(out2$sig_swe, paste0(outdir,"/sig_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$M1, paste0(outdir,"/M1_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$base1, paste0(outdir,"/base1_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$b1, paste0(outdir,"/b1_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out2$mid1, paste0(outdir,"/mid1_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$M2, paste0(outdir,"/M2_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$base2, paste0(outdir,"/base2_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$b2, paste0(outdir,"/b2_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$mid2, paste0(outdir,"/mid2_out_chain2_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	print(paste("end output",i))		

	}

	if(rn==3){	
	stan_model3 = stan(paste0(modDir), 
					data =  list(Nobs=dim(dat.swe5[dat.swe5$gcID==i,])[1], swe=dat.swe5$swe[dat.swe5$gcID==i], 
				day=(dat.swe5$jday[dat.swe5$gcID==i]-32)/(182-32),tempC=dat.swe5$tempCent[dat.swe5$gcID==i],vcf=dat.swe5$vcf[dat.swe5$gcID==i]),init=inits3,
				,chains=1, iter=3000)	
	print(paste("end model run",i))
	out3<- extract(stan_model3)	
	print(paste("extract variables",i))	
	write.table(out3$M0, paste0(outdir,"/M0_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$base0, paste0(outdir,"/base0_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$b0, paste0(outdir,"/b0_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$mid0, paste0(outdir,"/mid0_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	write.table(out3$sig_swe, paste0(outdir,"/sig_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$M1, paste0(outdir,"/M1_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$base1, paste0(outdir,"/base1_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$b1, paste0(outdir,"/b1_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out3$mid1, paste0(outdir,"/mid1_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	write.table(out1$M2, paste0(outdir,"/M2_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$base2, paste0(outdir,"/base2_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$b2, paste0(outdir,"/b2_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")
	write.table(out1$mid2, paste0(outdir,"/mid2_out_chain3_gc_",IDSglc$gcID[i],".csv"), sep=",")	
	print(paste("end output",i))		
				
	}
}

##this script ran on 806840 observations in the 6 classes in less than 2 hours