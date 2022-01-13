model{
  #folder run 4
  for(i in 1:Nobs){
    ####################
    #####likelihood#####
    ####################
    #slope of curve
    b0[i] ~ dnorm(mu.b0[i],tau.b0[glcYearID[i]])
    rep.b0[i] ~ dnorm(mu.b0[i],tau.b0[glcYearID[i]])
    #empirical regression
    mu.b0[i] <- betaB0[glcYearID[i]] + betaB1[glcYearID[i]]*TempAB[i] + betaB2[glcYearID[i]]*(CanopyB[i]-20) +
      betaB3[glcYearID[i]]*(sweDay[i]-107) + betaB4[glcYearID[i]]*(SweMax[i]-log(0.15))
    
  }
  ####################
  #####priors    #####
  ####################

  #### regression ####
  for(i in 1:NglcYear){
    #slope regression priors
    betaB0[i] ~ dnorm(mu.betaB0[glcID[i]],tau.betaB0[glcID[i]])
    betaB1[i] ~ dnorm(mu.betaB1[glcID[i]],tau.betaB1[glcID[i]])
    betaB2[i] ~ dnorm(mu.betaB2[glcID[i]],tau.betaB2[glcID[i]])
    betaB3[i] ~ dnorm(mu.betaB3[glcID[i]],tau.betaB3[glcID[i]])
    betaB4[i] ~ dnorm(mu.betaB4[glcID[i]],tau.betaB4[glcID[i]])
    
    #likelihood standard deviation
    tau.b0[i] <- pow(sig.b0[glcID[i]],-2)
    sig.b0[i] ~ dunif(0,100)
    
    
  }
  
  
  
  ####################
  ####hyper-priors####
  ####################
  #means
  for(i in 1:Nglc){
  mu.betaB0[i] ~ dnorm(0,0.00001)
  mu.betaB1[i] ~ dnorm(0,0.00001)
  mu.betaB2[i] ~ dnorm(0,0.00001)
  mu.betaB3[i] ~ dnorm(0,0.00001)
  mu.betaB4[i] ~ dnorm(0,0.00001)
  
  tau.betaB0[i] <- pow(sig.B0[i], -2)
  tau.betaB1[i] <- pow(sig.B1[i], -2)
  tau.betaB2[i] <- pow(sig.B2[i], -2)
  tau.betaB3[i] <- pow(sig.B3[i], -2)
  tau.betaB4[i] <- pow(sig.B4[i], -2)
  
  
  sig.B0[i] ~ dunif(0,1000)
  sig.B1[i] ~ dunif(0,1000)
  sig.B2[i] ~ dunif(0,1000)
  sig.B3[i] ~ dunif(0,1000)
  sig.B4[i] ~ dunif(0,1000)
  }



  
  ####################
  #####regression#####
  #####mean for  #####
  #####plotting  #####
  ####################
  for(i in 1:NglcYear){
    for(j in 1:200){
      mu.Temp[j,i] <- betaB0[i] + betaB1[i]*TempMean[j] 
      mu.Canopy[j,i] <- betaB0[i] + betaB2[i]*(CanopyMean[j]-20)
      mu.Onset[j,i] <- betaB0[i] + betaB3[i]*(SdayMean[j]-107)
      mu.Max[j,i] <- betaB0[i] + betaB4[i]*(MaxMean[j]-log(0.15))
    }
  }	
  

}
