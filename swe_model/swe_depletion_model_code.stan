data {
	int<lower=1> Nobs;
	real swe[Nobs];
	real day[Nobs];
	int<lower=1> Npixel;
	int<lower=1>  pixID[Nobs];
}
parameters{
		
	real<lower=0,upper=200> b0[Npixel]; 
	real<lower=0,upper=1> mid0[Npixel];	
	real<lower=0,upper=1> sig_swe;
	real<lower=0,upper=200> mu_b0; 
	real<lower=0,upper=1> mu_mid;
	real<lower=0,upper=200> sig_b0; 
	real<lower=0,upper=1> sig_mid;	
}	
model{


		
		sig_swe ~ uniform(0,1);
		mu_b0 ~ uniform(0,200);
		mu_mid ~ uniform(0,1);
		sig_b0 ~ uniform(0,200);
		sig_mid ~ uniform(0,1);
	
		b0 ~ normal(mu_b0,sig_b0)T(0,200);
		mid0 ~ normal(mu_mid,sig_mid)T(0,1);
		
	for(i in 1:Nobs){
	swe[i]~normal(1/(1+exp(b0[pixID[i]]*(day[i]-mid0[pixID[i]]))), sig_swe);
	
	}
}