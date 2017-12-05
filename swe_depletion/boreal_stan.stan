data {
	int<lower=1> Nobs;
	real swe[Nobs];
	int<lower=1> vegeC[Nobs];
	real day[Nobs];
	int<lower=1> Nveg;
}

parameters{
	real<lower=0> M[Nveg]; 
	real<lower=0> base[Nveg];		
	real<lower=0> b[Nveg]; 
	real<lower=0> sig_swe[Nveg];
	real<lower=1,upper= 1000> mu_M;
	real<lower=0,upper= 300> mu_base;
	real<lower=0,upper= 100> mu_b;
	real<lower=0,upper= 1000> sig_M;
	real<lower=0,upper= 100> sig_base;
	real<lower=0,upper= 100> sig_b;
}	
model{
	mu_M ~uniform(0,1000);
	mu_base ~uniform(0,300);
	mu_b ~uniform(0,100);
	sig_M ~ uniform(0,1000);
	sig_base~uniform(0,100);
	sig_b~uniform(0,100);

	for(i in 1:Nveg){
		M[i]~normal(mu_M, sig_M);
		base[i] ~normal(mu_base, sig_base);
		b[i] ~ normal(mu_b, sig_b);
		sig_swe[i] ~ uniform(0,1000);
	}	
	
	
	for(i in 1:Nobs){
	swe[i]~normal((M[vegeC[i]]/(1+exp(b[vegeC[i]]*(day[i]))))+base[vegeC[i]], sig_swe[vegeC[i]]);
	}



	
}