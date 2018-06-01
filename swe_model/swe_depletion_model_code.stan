data {
	int<lower=1> Nobs;
	real swe[Nobs];
	real day[Nobs];
	real tempC[Nobs];
	real vcf[Nobs];
}
parameters{
	real<lower=0,upper=1000> M0; 
	real<lower=0,upper=300> base0;		
	real<lower=0,upper=100> b0; 
	real<lower=0,upper=1> mid0;	
	real<lower=-1,upper=1> M1; 
	real<lower=-1,upper=1> base1;		
	real<lower=-1,upper=1> b1; 
	real<lower=-1,upper=1> mid1;		
	real<lower=0,upper=1000> sig_swe;
	real<lower=-1,upper=1> M2; 
	real<lower=-1,upper=1> base2;		
	real<lower=-1,upper=1> b2; 
	real<lower=-1,upper=1> mid2;		
}	
model{

		M0~uniform(0,1000);
		base0 ~uniform(0,300);
		b0 ~ uniform(0,100);
		mid0 ~ uniform(0,1);
		M1 ~ normal(0,.001)T[-1,1];
		base1 ~ normal(0,.001)T[-1,1];
		b1 ~ normal(0,.001)T[-1,1];
		mid1 ~ normal(0,.001)T[-1,1];
		sig_swe ~ uniform(0,1000);
		M2 ~ normal(0,.001)T[-1,1];
		base2 ~ normal(0,.001)T[-1,1];
		b2 ~ normal(0,.001)T[-1,1];
		mid2 ~ normal(0,.001)T[-1,1];
	
		
	for(i in 1:Nobs){
	swe[i]~normal(((M0+(M1*tempC[i])+(M2*vcf[i]))/(1+exp(((b0+(b1*tempC[i])+(b2*vcf[i]))*(day[i]-(mid0+(mid1*tempC[i])+(mid2*vcf[i])))))))+(base0+(base1*tempC[i])+(base2*vcf[i])), sig_swe);
	
	}
}