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

#run output directory
outdir <- "/local/swe/stan_all1"

#linux =1 or windows =2
runOS <- 2
#linux data directory first option, windows second optioon
DDdir <- c("/mnt/g/projects/boreal_swe_depletion/data/swe_depletion_model_data.csv",
				"z:\\projects\\boreal_swe_depletion\\data\\swe_depletion_model_data.csv")


+depletion<- function(b,day,mid,M,base){
+	(M/(1+exp(b*(day-mid))))+base
+}


#######################################################
# read in and organize data                           #
#######################################################
#read in data files
if(runOS==1){
	dat.swe <- read.csv(paste0(DDdir[1]))
}else{

	dat.swe <- read.csv(paste0(DDdir[2]))
}
print("finish reading in data")




#create a table of identifiers
IDScell <-
unique(dat.swe$glc1)


#take out any areas with no snow extent
#over the period
sweSums <- rowSums(sweS3[,3:94], na.rm=TRUE)
sweS4 <- sweS3[sweSums!=0,]

#see how many cells there are in each glc
gridC <- aggregate(sweS4[,1], by=list(sweS4[,2]), FUN="length")
colnames(gridC) <- c("glc", "count")

swe.vector <-as.vector(t(data.matrix(sweS4[,3:94])))
#turn into a data frame
sweDF <- data.frame(glc=rep(sweS4[,2],each=92),
			doy=rep(seq(61,152),times=dim(sweS4)[1]),
			gridID= rep(sweS4[,1],each=92),
			swe=swe.vector)

#plot the mean by vege type to see
sweM <- aggregate(sweDF$swe, by=list(sweDF$doy,sweDF$glc),FUN="mean", na.rm=TRUE)
colnames(sweM) <- c("doy","glc","swe")			
#see how many observations are in
plot(sweM$doy[sweM$glc==2],sweM$swe[sweM$glc==2])
plot(sweM$doy,sweM$swe)

#now subset to only include landcover classes intrested in
datgcI <- dat.gl[,1:2]
colnames(datgcI) <- c("glc","name")
gcInfo <- join(datgcI, gridC, by=c("glc"), type="left")
#don't include glc 3,10,15,16,18,20,21,22
gcUse <- gcInfo[gcInfo$glc!=2&gcInfo$glc!=3&gcInfo$glc!=15&gcInfo$glc!=16&gcInfo$glc!=10&gcInfo$glc!=14&
				gcInfo$glc!=9&gcInfo$glc!=14&gcInfo$glc<18,]
gcUse$gcID <- seq(1,  dim(gcUse)[1])
#subset to include only relevant land cover classes
sweDF2 <- join(sweDF, gcUse, by="glc", type="inner")
				
#get unique gridcells
gridID <- unique(data.frame(gridID=sweDF2$gridID, gcID=sweDF2$gcID))	
gridID$MgID <- seq(1, dim(gridID)[1])
#join back into
sweDF3 <- join(sweDF2, gridID, by=c("gridID","gcID"), type="left")

#######################################################
# set up model run                                    #
#######################################################



#inits for all possible runs					

inits2<-list(list(M=c(70),base=c(5),b=c(1),
	sig_swe=c(30)))
inits3<-list(list(M=c(60),base=c(10),b=c(.5)
		,sig_swe=c(10)))		
inits1<-list(list(M=c(50),base=c(1),b=c(.1),
			sig_swe=c(20)))
			
			
for(i in 1:dim(gcUse)[1]){		
			
	print(paste("start model run", i))			
	if(rn==1){		
	stan_model1 = stan("/home/hkropp/github/boreal_lw/swe_depletion/boreal_stan_nh_split.stan", 
					data = list(Nobs=dim(sweDF3[sweDF3$gcID==i,])[1], swe=sweDF3$swe[sweDF3$gcID==i], 
				day=sweDF3$doy[sweDF3$gcID==i]-(61+((152-61)/2))),init=inits1,
				,chains=1, iter=3000)	
	print(paste("end model run",i))	
	out1<- extract(stan_model1)
	print(paste("extract variables",i))	
	write.table(out1$M, paste0(outdir,"/M_out_chain1_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out1$base, paste0(outdir,"/base_out_chain1_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out1$b, paste0(outdir,"/b_out_chain1_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out1$sig_swe, paste0(outdir,"/sig_swe_out_chain1_gc_",gcUse$gcID[i],".csv"), sep=",")
	print("end output")				
				
	}
	if(rn==2){			
	stan_model2 = stan("/home/hkropp/github/boreal_lw/swe_depletion/boreal_stan_nh_split.stan", 
					data =  list(Nobs=dim(sweDF3[sweDF3$gcID==i,])[1], swe=sweDF3$swe[sweDF3$gcID==i], 
				day=sweDF3$doy[sweDF3$gcID==i]-(61+((152-61)/2))),init=inits2,
				,chains=1, iter=3000)	
	print(paste("end model run",i))
	out2<- extract(stan_model2)	
	print(paste("extract variables",i))	
	write.table(out2$M, paste0(outdir,"/M_out_chain2_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out2$base, paste0(outdir,"/base_out_chain2_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out2$b, paste0(outdir,"/b_out_chain2_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out2$sig_swe, paste0(outdir,"/sig_swe_out_chain2_gc_",gcUse$gcID[i],".csv"), sep=",")
	print(paste("end output",i))		

	}

	if(rn==3){	
	stan_model3 = stan("/home/hkropp/github/boreal_lw/swe_depletion/boreal_stan_nh_split.stan", 
					data =  list(Nobs=dim(sweDF3[sweDF3$gcID==i,])[1], swe=sweDF3$swe[sweDF3$gcID==i], 
				day=sweDF3$doy[sweDF3$gcID==i]-(61+((152-61)/2))),init=inits3,
				,chains=1, iter=3000)	
	print(paste("end model run",i))
	out3<- extract(stan_model3)	
	print(paste("extract variables",i))	
	write.table(out3$M, paste0(outdir,"/M_out_chain3_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out3$base, paste0(outdir,"/base_out_chain3_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out3$b, paste0(outdir,"/b_out_chain3_gc_",gcUse$gcID[i],".csv"), sep=",")
	write.table(out3$sig_swe, paste0(outdir,"/sig_swe_out_chain3_gc_",gcUse$gcID[i],".csv"), sep=",")
	print(paste("end output",i))		
				
	}
}

##this script ran on 806840 observations in the 6 classes in less than 2 hours