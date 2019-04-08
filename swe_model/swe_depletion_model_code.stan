data {
	int<lower=1> Nobs;
	int<lower=1> Nrep;
	real swe[Nobs];
	real day[Nobs];
	int<lower=1> Npixel;
	int<lower=1>  pixID[Nobs];
	real Rday[Nrep];
	int<lower=1>  RpixID[Nrep];
}
parameters{
		
	real<lower=0,upper=200> b0[Npixel]; 
	real<lower=0,upper=1> mid0[Npixel];	
	real<lower=0,upper=1> sig_swe;
	real<lower=0,upper=200> muB0; 
	real<lower=0,upper=1> muMid;
	real<lower=0,upper=200> sigB0; 
	real<lower=0,upper=1> sigMid;
	real swerep[Nrep];
}	
model{


		
		sig_swe ~ uniform(0,1);
		muB0 ~ uniform(0,200);
		muMid ~ uniform(0,1);
		sigB0 ~ uniform(0,200);
		sigMid ~ uniform(0,1);
	for(j in 1:Npixel){
		b0[j] ~ normal(muB0,sigB0)T[0,200];
		mid0[j] ~ normal(muMid,sigMid)T[0,1];
	}	
	for(i in 1:Nobs){
		swe[i]~normal(1/(1+exp(b0[pixID[i]]*(day[i]-mid0[pixID[i]]))), sig_swe);
	
	}
	for(m in 1:Nrep){
		swerep[m]~normal(1/(1+exp(b0[RpixID[m]]*(Rday[m]-mid0[RpixID[m]]))), sig_swe);
	}
}