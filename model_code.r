#########################################################
####code for a Bayesian mixed effects regression ########
####for snow water equivalent and tree cover     ########
#########################################################

model{
	####likelihood model 
	for(i in 1:Nobs){
		SWE[i]~dnorm(mu.SWE[i],tau.SWE)
		mu.SW[i]<-Beta1[LandID[i]]+Beta2[LandID[i]]*Tree.cov[i]+eps[yearID[i]]
	
	}
	#calculate identifiable intercept
	for(i in 1:Nland){
		Beta1star[i]<-Beta1[i]+eps.mean
	
	}

	####random effects 
	eps[1:Nyear]~dmnorm(mu.eps[1:Nyear], Omega.eps[,])
	
	for(i in 1:Nyear){
		mu.eps[i]<-0
		#calculate identifiable eps
		eps.star[i]<-eps[i]-eps.mean
		
	}
	eps.mean<-mean(eps[])
	
	
	#now calculate standard deviation for covariance matrix
	### Temporal covariance structure for year
	Omega.eps[1:NyearS,1:NyearS]<-inverse(Sigma.eps[1:NyearS,1:NyearS])
	for(y in 1:NyearS){
		for(m in 1:NyearS){
			Sigma.eps[y,m]<-(1/tau.eps)*exp(phi.eps*D.Y[y,m])
		
		}
	}
	#calculate the distance component of Sigma
	for(y in 1:NyearS){
		for(m in 1:NyearS){
			D.Y[y,m]<-sqrt(pow(xS[y]-xS[m],2)+ pow(yS[y]-yS[m],2))
		}
	}
	
	#set up priors for year covariance
	#use a folded t
	tau.eps<-pow(sig.eps,-2)
	sig.eps<-abs(t.eps)
	t.eps~dt(0,B,2)
	B<-1/(A*A)
	A<-10
	#prior for measure of autocorrelation
	phi.eps<-log(rho.eps)
	rho.eps~dbeta(1,1)
	
	####likelihood prior
	for(i in 1:Nland){
		Beta1[i]~dnorm(0,.001)
		Beta2[i]~dnorm(0,.001)
	}
	tau.SWE<-pow(sig.SWE, -2)
	sig.SWE~dunif(0, 1000)


}