#hierarchical regression of swe with covariates
#that also accounts for uncertainty in swe curve parameters
model{
	for(i in 1:Nobs){
	####################
	#####likelihood#####
	####################
		#slope of curve
		b0[i] ~ dnorm(mu.b0[i],tau.b0[glcIDB[i]])
		rep.b0[i] ~ dnorm(mu.b0[i],tau.b0[glcIDB[i]])
		#empirical regression
		mu.b0[i] <- betaB0[glcIDB[i]] + betaB1[glcIDB[i]]*TempAB[i] + betaB2[glcIDB[i]]*(CanopyB[i]-20) +
						 betaB3[glcIDB[i]]*(sweDay[i]-107) + betaB4[glcIDB[i]]*(Lat[i]-)
						 eps.b[GCyearB[i]] + eps.s[cellID[i]]
		#posterior predictive loss
		Sqdiff[i] <- pow(rep.b0[i] - b0[i],2)
		#log likelihood for waic
		loglike[i]<-.5*log(tau.b0[glcIDB[i]]/(2*3.141593))-((tau.b0[glcIDB[i]]/2)*pow(b0[i]-mu.b0[i],2))

	}
	####################
	#####priors    #####
	####################
	#### year random effects ####
	#year random effects
	for(i in 1:Ngcyear){
		#year random effects
		#with sweeping

		eps.b[i] ~ dnorm(0,tau.eb[ygcIDB[i]])

		eps.bS[i] <- eps.b[i] - epsb.bar[ygcIDB[i]]
	}
	for(i in 1:Nglc){
		#calculate means
		epsb.bar[i] <- mean(eps.b[startb[i]:endb[i]]) 	
		
		#variance terms for random effects
		tau.eb[i] <- pow(sig.eb[i],-2)
		sig.eb[i] ~ dgamma(0.0001,0.0001)
		
	}
	#### regression ####
	for(i in 1:Nglc){
		#slope regression priors
		betaB0[i] ~ dnorm(mu.betaB0,tau.betaB0)
		betaB1[i] ~ dnorm(mu.betaB1,tau.betaB1)
		betaB2[i] ~ dnorm(mu.betaB2,tau.betaB2)
		betaB3[i] ~ dnorm(mu.betaB3,tau.betaB3)

		
		#likelihood standard deviation
		tau.b0[i] <- pow(sig.b0[i],-2)
		sig.b0[i] ~ dunif(0,100)
		
		#calculate identifiable intercepts
		betaB0S[i] <- betaB0[i] + epsb.bar[i] + epsS.bar
		#look at intercept not on log scale
		trB0[i] <- exp(betaB0S[i])

	}

	
	####################
	#####cell      #####
	#####random    #####
	#####effects   #####
	####################

	for(j in 1:Ncell){
		#specify means
		eps.s[j] ~ dnorm(0,tau.es)
		#specify identifiable parameters
		eps.sS[j] <- eps.s[j]-epsS.bar
	}
		#variance terms for random effects
		tau.es <- pow(sig.es,-2)
		sig.es ~ dgamma(0.0001,0.0001)
		
	epsS.bar <- mean(eps.s[])	
	####################
	####hyper-priors####
	####################
	#means
	mu.betaB0 ~ dnorm(0,0.00001)
	mu.betaB1 ~ dnorm(0,0.00001)
	mu.betaB2 ~ dnorm(0,0.00001)
	mu.betaB3 ~ dnorm(0,0.00001)

	tau.betaB0 <- pow(sig.B0, -2)
	tau.betaB1 <- pow(sig.B1, -2)
	tau.betaB2 <- pow(sig.B2, -2)
	tau.betaB3 <- pow(sig.B3, -2)


	sig.B0 ~ dunif(0,1000)
	sig.B1 ~ dunif(0,1000)
	sig.B2 ~ dunif(0,1000)
	sig.B3 ~ dunif(0,1000)
	
	
	#Posterior predictive loss is the posterior mean of Dsum, must monitor Dsum
  Dsum <- sum(Sqdiff[])
	####################
	#####regression#####
	#####mean for  #####
	#####plotting  #####
	####################
	for(i in 1:Nglc){
		for(j in 1:200){
		mu.Temp[j,i] <- betaB0[i] + betaB1[i]*TempMean[j] 
		mu.Canopy[j,i] <- betaB0[i] + betaB2[i]*(CanopyMean[j]-20)
		mu.Onset[j,i] <- betaB0[i] + betaB3[i]*(SdayMean[j]-107)			 
		}
	}	

}
