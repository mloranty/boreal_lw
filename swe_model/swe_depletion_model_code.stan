data {
	int<lower=1> Nobs;
	real swe[Nobs];
	real day[Nobs];
}
parameters{
	real<lower=0,upper=1000> M; 
	real<lower=0,upper=300> base;		
	real<lower=0,upper=100> b; 
	real<lower=0,upper=1000> sig_swe;
	real<lower=0,upper=1> mid;
}	
model{

		M~uniform(0,1000);
		base ~uniform(0,300);
		b ~ uniform(0,100);
		sig_swe ~ uniform(0,1000);
		mid ~ uniform(0,1);
		
	for(i in 1:Nobs){
	swe[i]~normal((M/(1+exp(b*(day[i]-mid))))+base, sig_swe);
	}
}