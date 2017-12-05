library(coda)
library(plyr)
library(rstan)
outdir <- "/local/swe/stan_test1"

dat.swe <-read.csv("/local/swe/test_sub_swe.csv")

dat.gl <-read.csv("/local/swe/test_sub_glc.csv")

#######################################################
# set up model run                                    #
#######################################################
datalist <- list(Nobs=dim(dat.swe)[1], swe=dat.swe$swe, vegeC=dat.swe$gcID,
			day=dat.swe$doy,  mid=61+((152-61)/2),
			Nveg=dim(dat.gl)[1])
			
					

inits2<-list(list(M=c(70,70,70,70,70,70,70,70,70,70),base=c(5,5,5,5,5,5,5,5,5,5),b=c(1,1,1,1,1,1,1,1,1,1),
	sig.swe=c(30,30,30,30,30,30,30,30,30,30),mu.M=70,mu.base=5,mu.b=1,sig.M=20,sig.base=5,sig.b=.5))
inits3<-list(list(M=c(60,60,60,60,60,60,60,60,60,60),base=c(10,10,10,10,10,10,10,10,10,10),b=c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)
		,sig.swe=c(10,10,10,10,10,10,10,10,10,10),mu.M=60,mu.base=10,mu.b=.5,sig.M=30,sig.base=4,sig.b=.3))		
inits1<-list(list(M=c(50,50,50,50,50,50,50,50,50,50),base=c(1,1,1,1,1,1,1,1,1,1),b=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1),
			sig.swe=c(20,20,20,20,20,20,20,20,20,20),mu.M=50,mu.base=1,mu.b=.1,sig.M=10,sig.base=2,sig.b=.1))	
		
stan_modelh1 = stan("/home/hkropp/github/boreal_lw/swe_depletion/boreal_stan.stan", 
				data = list(Nobs=dim(dat.swe)[1], swe=dat.swe$swe, vegeC=dat.swe$gcID,
			day=dat.swe$doy-(61+((152-61)/2)),
			Nveg=dim(dat.gl)[1]),init=inits1,
			,chains=1, iter=3000)	

out1<- extract(stan_modelh1)			
			
			
stan_model2 = stan("/home/hkropp/github/boreal_lw/swe_depletion/boreal_stan_nh.stan", 
				data = list(Nobs=dim(dat.swe)[1], swe=dat.swe$swe, vegeC=dat.swe$gcID,
			day=dat.swe$doy-(61+((152-61)/2)),
			Nveg=dim(dat.gl)[1]),init=inits2,
			,chains=1, iter=3000)	

out2<- extract(stan_model2)	

stan_model3 = stan("/home/hkropp/github/boreal_lw/swe_depletion/boreal_stan_nh.stan", 
				data = list(Nobs=dim(dat.swe)[1], swe=dat.swe$swe, vegeC=dat.swe$gcID,
			day=dat.swe$doy-(61+((152-61)/2)),
			Nveg=dim(dat.gl)[1]),init=inits3,
			,chains=1, iter=3000)	

out3<- extract(stan_model1)				

write.table(out1$M, paste0(outdir,"/M_out_chain1.csv"), sep=",")
write.table(out2$M, paste0(outdir,"/M_out_chain2.csv"), sep=",")
write.table(out3$M, paste0(outdir,"/M_out_chain3.csv"), sep=",")



#view data
Mc1 <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\stan_test1\\stan_test1\\M_out_chain1.csv")
Mc2 <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\stan_test1\\stan_test1\\M_out_chain2.csv")
Mc3 <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\stan_test1\\stan_test1\\M_out_chain3.csv")
colnames(Mc1)<- paste0("M[",seq(1:10),"]")
colnames(Mc2)<- paste0("M[",seq(1:10),"]")
colnames(Mc3)<- paste0("M[",seq(1:10),"]")
Mtall <-mcmc.list( mcmc(Mc1),mcmc(Mc2),mcmc(Mc3))

plot(Mtall, ask=TRUE)