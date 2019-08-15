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
						 betaB3[glcIDB[i]]*(sweDay[i]-107)+ eps.b[GCyearB[i]] + eps.s[cellID[i]]
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
		#calculate identifiable intercepts
		betaB0S[i] <- betaB0[i] + epsb.bar[i] + epsS.bar
	


		
	}

		# variance parameters slope
		sig.vB ~ dgamma(0.0001,0.0001)		
	
	####################
	#####spatial   #####
	#####random    #####
	#####effects   #####
	####################
	eps.s[1:Ncell] ~ dmnorm(mu.epsS[1:Ncell],OmegaS[1:Ncell,1:Ncell])
	for(j in 1:Ncell){
		#specify means
		mu.epsS[j] <- 0
		#specify identifiable parameters
		eps.sS[j] <- eps.s[j]-epsS.bar
	}
	epsS.bar <- mean(eps.s[])
	#spatial covariance model for tree random effect
	OmegaS[1:Ncell,1:Ncell] <- inverse(SigmaS[1:Ncell,1:Ncell])
	#standard deviation for spatial covariance
	for(m in 1:Ncell){
			for(j in 1:Ncell){
				SigmaS[m,j] <- (1/tauS)*exp(phiS*DistS[m,j])
				#put distance in km rather than m
				DistS[m,j] <- sqrt(pow(x[j]-x[m],2)+ pow(y[m] - y[j], 2))/1000	
			}
		}	
	#priors for spatial covariance
	#folded t for standard deviation
	
		tauS <- pow(sigS,-2)
		sigS ~ dunif(0,100)
		#abs(t.S)
		#t.S ~ dt(0,p.S, 2)
		#p.S <- 1/(v.S*v.S)
		#v.S ~ dunif(0,100)
		
	#prior for autocorrelation
		phiS <-  log(rhoS)
		rhoS ~ dunif(0,1)
		#dbeta(alphaS,betaS)
		#alphaS ~ dunif(0,100)
		#betaS ~ dunif(0,100)
		
		
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


}