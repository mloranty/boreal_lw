data {
	int<lower=1> Nobs;
	real swe[Nobs];
	real day[Nobs];
	int<lower=1> Nveg;
}
parameters{
	real<lower=0,upper=1000> M; 
	real<lower=0,upper=300> base;		
	real<lower=0,upper=100> b; 
	real<lower=0,upper=1000> sig_swe;
}	
model{

		M~uniform(0,1000);
		base[i] ~uniform(0,300);
		b[i] ~ uniform(0,100);
		sig_swe[i] ~ uniform(0,1000);
		
	for(i in 1:Nobs){
	swe[i]~normal((M/(1+exp(b*(day[i]))))+base sig_swe);
	}
}