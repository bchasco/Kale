data{
  for(i in 1:(num_strata-1)){
    for(j in (i+1):num_strata){
      Rmat[i,j]<-0
      Rmat[j,i]<-0
    }
  }
  for(i in 1:num_strata){
    Rmat[i,i] <- 1
  }
}
model{
  #-----------------------
  #Derived parameters
  #-----------------------
  Sigma_p<-inverse(inv_Sigma_p)
  Sigma_pent<-inverse(inv_Sigma_pent)
  Sigma_phi<-inverse(inv_Sigma_phi)
  for(k in 1:num_strata){
    pent[1,k]<-1/(1+sum(delta[1:(s-1),k]))#similar to dirichlet trick except used additive log ratios
    psi[1,k]<-pent[1,k]
    lambda[s,k] <- 0
    psiPtot[k] <- sum(temp[1:s,k])
    Bstar[1:s,k] ~ dmulti(pent[1:s,k], Nsuper[k])
    for (i in 1:(s-1)){
      phi[i,k] <- ilogit(logit_phi[i,k]) ^ days[i+1,k]#days i+1 because delta 1 corresponds to period 2 (period 1 delta = 1)
      delta[i,k]<-exp(log_delta[i,k]) * days[i+1,k]#days i+1 because delta 1 corresponds to period 2 (period 1 delta = 1)
      psi[i+1,k] <- psi[i,k]*(1-p[i,k])*phi[i,k] + pent[i+1,k] *(phi[i,k]-1)/log(phi[i,k]) 
      lambda[i,k] <- phi[i,k]*(p[i+1,k]+(1-p[i+1,k])*lambda[i+1,k])
    }
    for(i in 2:s){
      pent[i,k]<-delta[i-1,k]/(1+sum(delta[1:(s-1),k])) #similar to dirichlet trick except used additive log ratios
    }
    for (i in 1:s){
      p[i,k] <- ilogit(logit_p[i,k])
      temp[i,k] <- psi[i,k]*p[i,k]
      multP[i,k] <- temp[i,k]/sum(temp[1:s,k])
      tau[i,k] <-p[i,k]/(p[i,k]+(1-p[i,k])*lambda[i,k]) 
    }
    #calculate process error variance and correlation matrix
    sigma_phi_process[k] <- sqrt(Sigma_phi[k,k])
    sigma_p_process[k] <- sqrt(Sigma_p[k,k])
    sigma_pent_process[k] <- sqrt(Sigma_pent[k,k])
    for (j in 1:num_strata){
      rho_phi[k,j] <- (Sigma_phi[k,j]/(sigma_phi_process[k]*sigma_phi_process[j]))
      rho_p[k,j] <- (Sigma_p[k,j]/(sigma_p_process[k]*sigma_p_process[j]))
      rho_pent[k,j] <- (Sigma_pent[k,j]/(sigma_pent_process[k]*sigma_pent_process[j]))
    }
  }
  #------
  #priors
  #------
  for(k in 1:num_strata){
    sigma_lambda[k] ~ dt(0,5^-2,1) T(0,)
    #nuisance variable to detect if n>0
    for (i in 1:s){
      v[i,k] ~ dbeta(a,b)
    }
    #Priors for first states in process model of prob capture, survival, birth:
    logit_phi[1,k] ~ dnorm(logit_phi_1_mu_prior, 5^-2)
    logit_p[1,k] ~ dnorm(0, 5^-2)
    log_delta[1,k] ~ dnorm(0,3^-2)
    #priors on abundance
    #priors on abundance
    Ntot[k] ~ dgamma(sigma_lambda[k]^-2,sigma_lambda[k]^-2)
    Nsuper[k] ~ dpois(Ntot[k])
  }
  #process model priors
  for (i in 2:(s-1)){
    logit_phi[i,1:num_strata] ~ dmnorm(logit_phi[i-1,1:num_strata], inv_Sigma_phi)
    log_delta[i,1:num_strata] ~ dmnorm(log_delta[i-1,1:num_strata], inv_Sigma_pent)#similar to dirichlet trick except used additive log ratios
  }
  for(i in 2:s){
    logit_p[i,1:num_strata] ~ dmnorm(logit_p[i-1,1:num_strata], inv_Sigma_p)
  }
  #priors for process error covariance matrices
  inv_Sigma_p ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  inv_Sigma_phi ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  inv_Sigma_pent ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  #-----------
  #Likelihoods
  #-----------
  for(k in 1:num_strata){
    uTot[k] ~ dbin(psiPtot[k],Nsuper[k])
    u[1:s, k] ~ dmulti(multP[1:s,k],uTot[k])
    for (i in 1:(s-1)){
      R[i,k] ~ dbin(v[i,k],n[i,k])
      r[i,k] ~ dbin(lambda[i,k],R[i,k])
    }
    for (i in 2:(s-1)){
      m[i,k] ~dbin(tau[i,k],T[i,k])
    }
  }
}
