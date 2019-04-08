#hierarchical regression of swe with covariates
#that also accounts for uncertainty in swe curve parameters
model{
	for(i in 1:Nobs){
		#midpoint of curve
		maxD[i] ~ dnorm(mu.max[i],tau.max)
		rep.max[i] ~ dnorm(mu.max[i],tau.max)
		#empirical regression
		mu.max[i] <- betaM0[glcIDM[i]] + betaM1[glcIDM[i]]*TempAM[i] + betaM2[glcIDM[i]]*CanopyM[i] +
						 betaM3[glcIDM[i]]*(Lat[i]-55)+ eps.max[GCyearM[i]] 


		
		#slope of curve
		b0[i] ~ dnorm(mu.b0[i],tau.b0[i])
		rep.b0[i] ~ dnorm(mu.b0[i],tau.b0[i])
		#empirical regression
		mu.b0[i] <- betaB0[glcIDB[i]] + betaB1[glcIDB[i]]*TempAB[i] + betaB2[glcIDB[i]]*CanopyB[i] +
						 betaB3[glcIDB[i]]*(sweMaxB[i]-0.1)+ eps.b[GCyearB[i]] 
		#error model 
		tau.b0[i] <- pow(sig.b0[i],-2)
		sig.b0[i] <- sig.modB[i] + sig.vB


	}
	#priors
	for(i in 1:Ngcyear){
	#year random effects
	#with sweeping
		eps.max[i] ~ dnorm(0,tau.em[ygcIDM[i]])
		eps.b[i] ~ dnorm(0,tau.eb[ygcIDB[i]])
		eps.maxS[i] <- eps.max[i] - epsm.bar[ygcIDM[i]]
		eps.bS[i] <- eps.b[i] - epsb.bar[ygcIDB[i]]
	}
	for(i in 1:Nglc){
		#midpoint regression priors
		betaM0[i] ~ dnorm(mu.betaM0,tau.betaM0)
		betaM1[i] ~ dnorm(mu.betaM1,tau.betaM1)
		betaM2[i] ~ dnorm(mu.betaM2,tau.betaM2)
		betaM3[i] ~ dnorm(mu.betaM3,tau.betaM3)

		#slope regression priors
		betaB0[i] ~ dnorm(mu.betaB0,tau.betaB0)
		betaB1[i] ~ dnorm(mu.betaB1,tau.betaB1)
		betaB2[i] ~ dnorm(mu.betaB2,tau.betaB2)
		betaB3[i] ~ dnorm(mu.betaB3,tau.betaB3)
		#calculate identifiable intercepts
		betaB0S[i] <- betaB0[i] + epsb.bar[i]
		betaM0S[i] <- betaM0[i] + epsm.bar[i]
		#calculate means
		epsb.bar[i] <- mean(eps.b[startb[i]:endb[i]]) 	
		epsm.bar[i] <- mean(eps.max[startm[i]:endm[i]]) 
		#variance terms for random effects
		tau.em[i] <- pow(sig.em[i],-2)
		tau.eb[i] <- pow(sig.eb[i],-2)
		
		sig.em[i] ~ dgamma(0.0001,0.0001)
		sig.eb[i] ~ dgamma(0.0001,0.0001)
	}

	
		# variance parameters midpoint
		tau.max <- pow(sig.vM,-2)
		sig.vM ~ dgamma(0.0001,0.0001)
		# variance parameters slope
		sig.vB ~ dgamma(0.0001,0.0001)		
	
	#hyper priors
	#means
	mu.betaM0 ~ dnorm(0,0.00001)
	mu.betaM1 ~ dnorm(0,0.00001)
	mu.betaM2 ~ dnorm(0,0.00001)
	mu.betaM3 ~ dnorm(0,0.00001)


	mu.betaB0 ~ dnorm(0,0.00001)
	mu.betaB1 ~ dnorm(0,0.00001)
	mu.betaB2 ~ dnorm(0,0.00001)
	mu.betaB3 ~ dnorm(0,0.00001)

	
	#standard deviation
	tau.betaM0 <- pow(sig.M0, -2)
	tau.betaM1 <- pow(sig.M1, -2)
	tau.betaM2 <- pow(sig.M2, -2)
	tau.betaM3 <- pow(sig.M3, -2)

	
	sig.M0 ~ dunif(0,1000)
	sig.M1 ~ dunif(0,1000)
	sig.M2 ~ dunif(0,1000)
	sig.M3 ~ dunif(0,1000)



	tau.betaB0 <- pow(sig.B0, -2)
	tau.betaB1 <- pow(sig.B1, -2)
	tau.betaB2 <- pow(sig.B2, -2)
	tau.betaB3 <- pow(sig.B3, -2)

	sig.B0 ~ dunif(0,1000)
	sig.B1 ~ dunif(0,1000)
	sig.B2 ~ dunif(0,1000)
	sig.B3 ~ dunif(0,1000)

}