#simplest possible model
model{
	for(i in 1:Nobs){
	swe[i]~dnorm(mu.swe[i], tau.swe[vegeC[i]])
	mu.swe[i] <- (M[vegeC[i]]/(1+exp(b[vegeC[i]]*(day[i]-mid))))+base[vegeC[i]]
	}

	for(i in 1:NVeg){
		M[i]~dnorm(mu.M, tau.M)T(0,)
		base[i] ~dnorm(mu.base, tau.base)T(0,)
		b[i] ~ dnorm(mu.b, tau.b)T(0,)
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