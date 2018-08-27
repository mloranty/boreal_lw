data {
	int<lower=1> Nobs;	
	real swe[Nobs];
	real day[Nobs];
	real temp[Nobs];
	real treeCov[Nobs];
	

}
parameters{
		
	real<lower=-20,upper=5> beta0; 
	real beta1; 
	real beta2; 
	real beta3; 
	
	real<lower=-20,upper=0> alpha0; 
	real alpha1; 
	real alpha2; 
	real alpha3; 	
	real<lower=0,upper=1> sig_swe;

	
}	
model{


		beta0 ~ uniform(-20,5);
		beta1 ~ normal(0,10);
		beta2 ~ normal(0,10);
		beta3 ~ normal(0,10);
		
		alpha0 ~ uniform(-20,0);
		alpha1 ~ normal(0,1);
		alpha2 ~ normal(0,1);
		alpha3 ~ normal(0,1);		
	
		
		sig_swe ~ uniform(0,1);

	for(i in 1:Nobs){
	swe[i]~ normal(1/(1+exp((exp(beta0 + (beta1*temp[i]) + (beta2*treeCov[i]))) *(day[i]-(exp(alpha0 + (alpha1*temp[i]) + (alpha2*treeCov[i])))))), sig_swe);
	
	}

}

