data {
	int<lower=1> Nobs;
	real swe[Nobs];
	int<lower=1> vegeC[Nobs];
	real day[Nobs];
	int<lower=1> Nveg;
}
parameters{
	real<lower=0,upper=1000> M[Nveg]; 
	real<lower=0,upper=300> base[Nveg];		
	real<lower=0,upper=100> b[Nveg]; 
	real<lower=0,upper=1000> sig_swe[Nveg];
}	
model{
	for(i in 1:Nveg){
		M[i]~uniform(0,1000);
		base[i] ~uniform(0,300);
		b[i] ~ uniform(0,100);
		sig_swe[i] ~ uniform(0,1000);
	}	
	for(i in 1:Nobs){
	swe[i]~normal((M[vegeC[i]]/(1+exp(b[vegeC[i]]*(day[i]))))+base[vegeC[i]], sig_swe[vegeC[i]]);
	}
}