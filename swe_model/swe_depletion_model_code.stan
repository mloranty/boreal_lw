data {
	int<lower=1> Nobs;
	real swe[Nobs];
	real day[Nobs];
	real tempC[Nobs];
	real vcf[Nobs];
}
parameters{
		
	real<lower=0,upper=100> b0; 
	real<lower=0,upper=1> mid0;	
	
}	
model{

	for(j in 1:Npixel){
		b0[j] ~ uniform(0,100);
		mid0[j] ~ uniform(0,1);

	}	
	for(i in 1:Nobs){
	swe[i]~normal(1/(1+exp(b0[pixID[i]]*(day-mid0[pixID[i]])), sig_swe);
	
	}
}