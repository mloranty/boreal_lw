#hierarchical regression of swe with covariates
#that also accounts for uncertainty in swe curve parameters
model{
	for(i in 1:Nobs){
	####################
	#####likelihood#####
	####################
		#slope of curve
		b0[i] ~ dnorm(mu.b0[i],tau.b0[i])
		rep.b0[i] ~ dnorm(mu.b0[i],tau.b0[i])
		#empirical regression
		mu.b0[i] <- betaB0[glcIDB[i]] + betaB1[glcIDB[i]]*TempAB[i] + betaB2[glcIDB[i]]*CanopyB[i] +
						 betaB3[glcIDB[i]]*(sweDay[i]-107) + betaB4[glcIDB[i]]*TempAB[i]*CanopyB[i] +
						 betaB5[glcIDB[i]]*TempAB[i]*(sweDay[i]-107) + betaB6[glcIDB[i]]*CanopyB[i]*(sweDay[i]-107) +
						 eps.b[GCyearB[i]] + eps.s[cellID[i]]
		#error model 
		tau.b0[i] <- pow(sig.b0[i],-2)
		sig.b0[i] <- sig.modB[i] + sig.vB
		#posterior predictive loss
		Sqdiff[i] <- pow(rep.b0[i] - b0[i],2)
		#log likelihood for waic
		loglike[i]<-.5*log(tau.b0[i]/(2*3.141593))-((tau.b0[i]/2)*pow(b0[i]-mu.b0[i],2))

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
		betaB4[i] ~ dnorm(mu.betaB4,tau.betaB4)
		betaB5[i] ~ dnorm(mu.betaB5,tau.betaB5)
		betaB6[i] ~ dnorm(mu.betaB6,tau.betaB6)
		#calculate identifiable intercepts
		betaB0S[i] <- betaB0[i] + epsb.bar[i] + epsS.bar
		
		#estimate effects with interactions
		#low canopy is 10% and temp is -2 and melt doy is 60
		covEffectLow[i] <- betaB2[i] + (betaB4[i]*-2 ) + betaB6[i]*(60-107)
		TempEffectLow[i] <-  betaB1[i] + betaB4[i]*10 +  betaB5[i]*(60-107)
		MeltEffectLow[i] <-  betaB3[i] + (betaB5[i]*-2) + betaB6[i]*10
		#mid canopy is 30% and temp is 1 and melt doy is 100
		covEffectMid[i] <- betaB2[i] + (betaB4[i]*1 ) + betaB6[i]*(100-107)
		TempEffectMid[i] <-  betaB1[i] + betaB4[i]*30 +  betaB5[i]*(100-107)
		MeltEffectMid[i] <-  betaB3[i] + (betaB5[i]*1) + betaB6[i]*30
		#high canopy is 60% and temp is 5 and melt doy is 140

		covEffectHigh[i] <- betaB2[i] + (betaB4[i]*5 ) + betaB6[i]*(140-107)
		TempEffectHigh[i] <-  betaB1[i] + betaB4[i]*60 +  betaB5[i]*(140-107)
		MeltEffectHigh[i] <-  betaB3[i] + (betaB5[i]*5) + betaB6[i]*60
		
		#look at when one interaction is zero
		#and at low 
		#cover with low temp but middle of melt period
		covTempLow[i] <- betaB2[i] + (betaB4[i]*-2 )
		#cov with early melt but temp at zero
		covMeltLow[i] <- betaB2[i] + betaB6[i]*(60-107)
		#melt with a low temp and tree cover of zero
		meltTempLow[i] <- betaB3[i] + (betaB5[i]*-2) 
		#melt day with low canopy cover and temp of zero
		meltCovLow[i] <-  betaB3[i]  + betaB6[i]*10
		#temperature with low canopy cover
		tempCovLow[i] <- betaB1[i] + betaB4[i]*10 
		#temperature with early melt day
		tempMeltLow[i] <- betaB1[i] +  betaB5[i]*(60-107)
		#look at when one interaction is zero
		#and at high
		#cover with high temp but middle of melt period
		covTempHigh[i] <- betaB2[i] + (betaB4[i]*5 )
		#cov with late melt but temp at zero
		covMeltHigh[i] <- betaB2[i] + betaB6[i]*(140-107)
		#melt with a high temp and tree cover of zero
		meltTempHigh[i] <- betaB3[i] + (betaB5[i]*5) 
		#melt day with high canopy cover and temp of zero
		meltCovHigh[i] <-  betaB3[i]  + betaB6[i]*60
		#temperature with high canopy cover
		tempCovHigh[i] <- betaB1[i] + betaB4[i]*60 
		#temperature with late melt day and no cnaopy
		tempMeltHigh[i] <- betaB1[i] +  betaB5[i]*(140-107)
		
	}

		# variance parameters slope
		sig.vB ~ dgamma(0.0001,0.0001)		
	
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
	mu.betaB4 ~ dnorm(0,0.00001)
	mu.betaB5 ~ dnorm(0,0.00001)
	mu.betaB6 ~ dnorm(0,0.00001)

	tau.betaB0 <- pow(sig.B0, -2)
	tau.betaB1 <- pow(sig.B1, -2)
	tau.betaB2 <- pow(sig.B2, -2)
	tau.betaB3 <- pow(sig.B3, -2)
	tau.betaB4 <- pow(sig.B4, -2)
	tau.betaB5 <- pow(sig.B5, -2)
	tau.betaB6 <- pow(sig.B6, -2)

	sig.B0 ~ dunif(0,1000)
	sig.B1 ~ dunif(0,1000)
	sig.B2 ~ dunif(0,1000)
	sig.B3 ~ dunif(0,1000)
	sig.B4 ~ dunif(0,1000)
	sig.B5 ~ dunif(0,1000)
	sig.B6 ~ dunif(0,1000)
	
	
	#Posterior predictive loss is the posterior mean of Dsum, must monitor Dsum
  Dsum <- sum(Sqdiff[])


}