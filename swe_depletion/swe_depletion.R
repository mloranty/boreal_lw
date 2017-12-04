


dat.swe <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\glc_SWE_M2BGS_ease2_2010.csv")


depletion<- function(b,day,mid,M){
	(M/(1+exp(b*(day-mid))))+2
}

plot(seq(50,100),depletion(1,seq(50,100),75,10), type="l")
