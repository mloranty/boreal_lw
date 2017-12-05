library(coda)
library(plyr)
library(rstan)


#start with a test year
dat.swe <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\test_sub_swe.csv")
dat.gl <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\test_sub_glc.csv")

#######################################################
# set up model run                                    #
#######################################################
datalist <- list(Nobs=dim(dat.swe)[1], swe=dat.swe$swe, vegeC=dat.swe$gcID,
			day=dat.swe$doy,  mid=61+((152-61)/2),
			Nveg=dim(dat.gl)[1])
			
					
inits<-list(list(M=c(50,50,50,50,50,50,50,50,50,50),		
			base=c(1,1,1,1,1,1,1,1,1,1),
			b=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1),
			sig.swe=c(20,20,20,20,20,20,20,20,20,20),
			mu.M=50,mu.base=1,mu.b=.1,sig.M=10,sig.base=2,
			sig.b=.1),
			list(M=c(70,70,70,70,70,70,70,70,70,70),		
			base=c(5,5,5,5,5,5,5,5,5,5),
			b=c(1,1,1,1,1,1,1,1,1,1,1),
			sig.swe=c(30,30,30,30,30,30,30,30,30,30),
			mu.M=70,mu.base=5,mu.b=1,sig.M=20,sig.base=5,
			sig.b=.5),
			list(M=c(60,60,60,60,60,60,60,60,60,60),		
			base=c(10,10,10,10,10,10,10,10,10,10),
			b=c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5),
			sig.swe=c(10,10,10,10,10,10,10,10,10,10),
			mu.M=60,mu.base=10,mu.b=.5,sig.M=30,sig.base=4,
			sig.b=.3))		
			
stan_model = stan("c:\\Users\\hkropp\\Documents\\test_stan\\boreal\\boreal_stan.stan", 
				data = datalist,init=inits,
					chains=3, verbose=TRUE)			
			
			
parms <- c("M","base","sig.swe","b","sig.swe","sig.bV",
			"mu.base","mu.b",
			"sig.M","sig.base","sig.b")		