data {
	int<lower=1> Nobs;	
	real swe[Nobs];
	real day[Nobs];
	real temp[Nobs];
	real treeCov[Nobs];
	real year[Nobs];
	

}
parameters{
		
	real<lower=0,upper=200> beta0; 
	real<lower=-1,upper=1> beta1; 
	real<lower=-1,upper=1> beta2; 
	real<lower=-1,upper=1> beta3; 
	
	real<lower=0,upper=1> alpha0; 
	real<lower=-.1,upper=.1> alpha1; 
	real<lower=-.1,upper=.1> alpha2; 
	real<lower=-.1,upper=.1> alpha3; 	
	real<lower=0,upper=1> sig_swe;

	
}	
model{


		beta0 ~ uniform(0,100);
		beta1 ~ normal(0,1)T[-1,1];;
		beta2 ~ normal(0,1)T[-1,1];;
		beta3 ~ normal(0,1)T[-1,1];;
		
		alpha0 ~ uniform(0,1);
		alpha1 ~ normal(0,.1)T[-.1,.1];;
		alpha2 ~ normal(0,.1)T[-.1,.1];;
		alpha3 ~ normal(0,.1)T[-.1,.1];;		
	
		
		sig_swe ~ uniform(0,1);

	for(i in 1:Nobs){
	swe[i]~ normal(1/(1+exp((beta0 + (beta1*temp[i]) + (beta2*treeCov[i])) *(day[i]-(alpha0 + (alpha1*temp[i]) + (alpha2*treeCov[i]))))), sig_swe);
	
	}

}

