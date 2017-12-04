#simplest possible model
model{
	for(i in 1:Nobs){
	swe[i]~dnorm(mu.swe[i], tau.swe[vegeC[i]])
	mu.swe[i] <- (M[gridC[i]]/(1+exp(b[gridC[i]]*(day[i]-mid))))+base[gridC[i]]
	}

	for(i in 1:NGridC){
		M[i]~dnorm(mu.MV[vege[i]], tau.MV[vege[i]])T(0,)
		base[i] ~dnorm(mu.baseV[vege[i]], tau.baseV[vege[i]])T(0,)
		b[i] ~ dnorm(mu.bV[vege[i]], tau.bV[vege[i]])T(0,)
		
	}
	for(i in 1:NVeg){
		mu.MV[i] ~dnorm(mu.M, tau.M)T(0,)
		mu.baseV[i] ~dnorm(mu.base, tau.base)T(0,)
		mu.bV[i] ~ dnorm(mu.b, tau.b)T(0,)
		tau.MV[i] <- pow(sig.MV[i],-2)
		tau.baseV[i] <- pow(sig.baseV[i],-2)
		tau.bV[i] <- pow(sig.bV[i],-2)
		sig.MV[i] ~ dunif(0,1000)
		sig.baseV[i]~dunif(0,100)
		sig.bV[i]~dunif(0,100)
		tau.swe[i] <- pow(sig.swe[i],-2)
		sig.swe[i] ~ dunif(0,1000)
	}
	mu.M ~dunif(0,1000)
	mu.base ~dunif(0,300)
	mu.b ~dunif(0.00001,100)
	tau.M <- pow(sig.M,-2)
	tau.base <- pow(sig.base,-2)
	tau.b <- pow(sig.b,-2)
	sig.M ~ dunif(0,1000)
	sig.base~dunif(0,100)
	sig.b~dunif(0,100)
}