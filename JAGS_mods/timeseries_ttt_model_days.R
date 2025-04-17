model{
  #-----------------------
  #process model equations
  #-----------------------
  pent[1]<-1/(1+sum(delta[1:(s-1)]))#similar to dirichlet trick except used additive log ratios
  psi[1]<-pent[1]
  lambda[s] <- 0
  psiPtot <- sum(temp[1:s])
  Bstar[1:s] ~ dmulti(pent[1:s], Nsuper)
  for (i in 1:(s-1)){
    phi[i] <- ilogit(logit_phi[i]) ^ days[i+1]
    delta[i] <- exp(log_delta[i]) * days[i+1]
    psi[i+1] <- psi[i]*(1-p[i])*phi[i] + pent[i+1] *(phi[i]-1)/log(phi[i]) 
    lambda[i] <- phi[i]*(p[i+1]+(1-p[i+1])*lambda[i+1])
  }
  for(i in 2:s){
    pent[i]<-delta[i-1]/(1+sum(delta[1:(s-1)])) #similar to dirichlet trick except used additive log ratios
  }
  for (i in 1:s){
    p[i] <- ilogit(logit_p[i])
    temp[i] <- psi[i]*p[i]
    multP[i] <- temp[i]/sum(temp[1:s])
    tau[i] <-p[i]/(p[i]+(1-p[i])*lambda[i]) 
  }
  #------
  #priors
  #------
  #nuisance variable to detect if R>0
  for (i in 1:s){
    v[i] ~ dbeta(a,b)
  }
  #process model priors
  for (i in 2:(s-1)){
    logit_phi[i] ~ dnorm(logit_phi[i-1], sigma_phi^-2)
    log_delta[i] ~ dnorm(log_delta[i-1], sigma_delta^-2)#similar to dirichlet trick except used additive log ratios
  }
  for(i in 2:s){
    logit_p[i] ~ dnorm(logit_p[i-1], sigma_p^-2)
  }
  #Priors for first states in process model of prob birth, capture, survival:
  logit_phi[1] ~ dnorm(logit_phi_1_mu_prior,5^-2)
  logit_p[1] ~ dnorm(0, 5^-2)
  log_delta[1] ~ dnorm(log(1), 5^-2)#similar to dirichlet trick except used additive log ratios
  #priors for process error sd
  sigma_lambda ~ dt(0,2.5^-2,1) T(0,)
  sigma_phi ~ dt(0,2.5^-2,1) T(0,)
  sigma_p ~ dt(0,2.5^-2,1) T(0,)
  sigma_delta ~ dt(0,2.5^-2,1) T(0,)
  #priors on abundance
  Ntot ~ dgamma(sigma_lambda^-2,sigma_lambda^-2)
  Nsuper ~ dpois(Ntot)
  #-----------
  #Likelihoods
  #-----------
  uTot ~ dbin(psiPtot,Nsuper)
  u[1:s] ~ dmulti(multP[1:s],uTot)
  for (i in 1:(s-1)){
    R[i] ~ dbin(v[i],n[i])
    r[i] ~ dbin(lambda[i],R[i])
  }
  for (i in 2:(s-1)){
    m[i] ~dbin(tau[i],T[i])
  }
}
