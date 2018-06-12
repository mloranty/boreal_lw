data {
	int<lower=1> Nobs;
	real swe[Nobs];
	real day[Nobs];
	int<lower=1> Npixel;
	int<lower=1>  pixID[Nobs];
}
parameters{
		
	real<lower=0,upper=100> b0[Npixel]; 
	real<lower=0,upper=1> mid0[Npixel];	
	real<lower=0,upper=1> sig_swe;
}	
model{


		b0 ~ uniform(0,100);
		mid0 ~ uniform(0,1);
		sig_swe ~ uniform(0,1);

	for(i in 1:Nobs){
	swe[i]~normal(1/(1+exp(b0[pixID[i]]*(day[i]-mid0[pixID[i]]))), sig_swe);
	
	}
}