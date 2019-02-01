#hierarchical regression of swe with covariates
#that also accounts for uncertainty in swe curve parameters
model{
	for(i in 1:Nobs){
		#midpoint of curve
		mid[i] ~ dnorm(mu.mid[i],tau.mid[glcIDM[i]])
		#empirical regression
		mu.mid[i] <- betaM0[glcIDM[i]] + betaM1[glcIDM[i]]*TempAM[i] + betaM2[glcIDM[i]]*CanopyM[i] +
						 betaM3[glcIDM[i]]*(yearM[i]-2000)
		#error model and standard deviation
		tau.mid[i] <- tau.modM[i] + tau.vM[glcIDM[i]]
		tau.modM[i] <- pow(sig.modM[i],-2)
		
		#slope of curve
		b0[i] ~ dnorm(mu.b0[i],tau.b0[glcIDM[i]])
		#empirical regression
		mu.b0[i] <- betaB0[glcIDB[i]] + betaB1[glcIDB[i]]*TempAB[i] + betaB2[glcIDB[i]]*CanopyB[i] +
						 betaB3[glcIDB[i]]*(yearB[i]-2000)
		#error model 
		tau.b0[i] <- tau.modB[i] + tau.vB[glcIDB[i]]
		tau.modB[i] <- pow(sig.modB[i],-2)

	}
	#priors
	for(i in 1:Nglc){
		#midpoint regression priors
		betaM0[i] ~ dnorm(0,0.0000001)
		betaM1[i] ~ dnorm(0,0.0000001)
		betaM2[i] ~ dnorm(0,0.0000001)
		betaM3[i] ~ dnorm(0,0.0000001)

		#slope regression priors
		betaB0[i] ~ dnorm(0,0.0000001)
		betaB1[i] ~ dnorm(0,0.0000001)
		betaB2[i] ~ dnorm(0,0.0000001)
		betaB3[i] ~ dnorm(0,0.0000001)

		# variance parameters midpoint
		tau.vM[i] <- pow(sig.vM[i],-2)
		sig.vM[i] ~ dunif(0,10)
		# variance parameters slope
		tau.vB[i] <- pow(sig.vB[i],-2)
		sig.vB[i] ~ dunif(0,500)		
	}


}