model{
###################################################################################################################### #
###########              Multivariate, time-series JS (ttt) model                                           ###########
####################################################################################################################### #   
  #-----------------------
  #Derived parameters
  #-----------------------
  Sigma_p<-inverse(inv_Sigma_p)
  Sigma_phi<-inverse(inv_Sigma_phi)
  Sigma_pent<-inverse(inv_Sigma_pent)
  for(k in 1:num_strata){
    #pent[1,k]<-1/(1+sum(delta[1:(s-1),k]))#similar to dirichlet trick except used additive log ratios #TB4.13.2020
    psi[1,k]<-pent[1,k]
    lambda[s,k] <- 0
    psiPtot[k] <- sum(temp[1:s,k])
    Bstar[1:s,k] ~ dmulti(pent[1:s,k], Nsuper[k])
    for (i in 1:(s-1)){
      phi[i,k] <- ilogit(logit_phi[i,k]) ^ time[i]
      #delta[i,k]<-exp(log_delta[i,k]) * time[i] #TB4.13.2020
      psi[i+1,k] <- psi[i,k]*(1-p[i,k])*phi[i,k] + pent[i+1,k] *(phi[i,k]-1)/log(phi[i,k]) 
      lambda[i,k] <- phi[i,k]*(p[i+1,k]+(1-p[i+1,k])*lambda[i+1,k])
    }
    #for(i in 2:s){ #TB4.13.2020
      #pent[i,k]<-delta[i-1,k]/(1+sum(delta[1:(s-1),k])) #similar to dirichlet trick except used additive log ratios
    #}
    for (i in 1:s){
      delta[i,k]<-exp(log_delta[i,k]) * time[i] #TB4.13.2020
      pent[i,k]<-delta[i,k]/sum(delta[1:s,k]) #similar to dirichlet trick except used additive log ratios  #TB4.13.2020
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
    sigma_lambda[k] ~ dt(hyper_value_sigma_lambda_mu, hyper_value_sigma_lambda_sd^-2, 1) T(0,)
    #nuisance variable to detect if n>0
    for (i in 1:s){
      v[i,k] ~ dbeta(hyper_value_beta_v, hyper_value_beta_v)
    }
    #Priors for first states in process model of prob capture, survival, birth:
    logit_phi[1,k] ~ dnorm(hyper_value_logit_phi_1_mu, hyper_value_logit_phi_1_sd^-2)
    logit_p[1,k]   ~ dnorm(hyper_value_logit_p_1_mu,   hyper_value_logit_p_1_sd^-2)
    log_delta[1,k] ~ dnorm(hyper_value_log_delta_1_mu, hyper_value_log_delta_1_sd^-2)
    #priors on abundance
    Ntot[k] ~ dgamma(sigma_lambda[k]^-2, sigma_lambda[k]^-2)
    Nsuper[k] ~ dpois(Ntot[k])
  }
  #process model priors
  for (i in 2:(s-1)){
    logit_phi[i,1:num_strata] ~ dmnorm(logit_phi[i-1,1:num_strata], inv_Sigma_phi)
    #log_delta[i,1:num_strata] ~ dmnorm(log_delta[i-1,1:num_strata], inv_Sigma_pent)#similar to dirichlet trick except used additive log ratios TB4.13.2020
  }
  for(i in 2:s){
    log_delta[i,1:num_strata] ~ dmnorm(log_delta[i-1,1:num_strata], inv_Sigma_pent)#similar to dirichlet trick except used additive log ratios TB4.13.2020
    logit_p[i,1:num_strata] ~ dmnorm(logit_p[i-1,1:num_strata], inv_Sigma_p)
  }
  #priors for process error covariance matrices
  inv_Sigma_p    ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  inv_Sigma_phi  ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  inv_Sigma_pent ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  #-----------
  #Likelihoods
  #-----------
  for(k in 1:num_strata){
    uTot[k] ~ dbin(psiPtot[k],Nsuper[k])
    u[1:s, k] ~ dmulti(multP[1:s,k],uTot[k])
    for (i in 1:(s-1)){
      R[i,k] ~ dbin(v[i,k], n[i,k])
      r[i,k] ~ dbin(lambda[i,k],R[i,k])
    }
    for (i in 2:(s-1)){
      m[i,k] ~dbin(tau[i,k],T[i,k])
    }
  }

####################################################################################################################### #
###########   AUC, RT, ESC from Manske & Schwarz 2000
####################################################################################################################### #
  for(k in 1:num_strata){
    B[1,k] <- Bstar[1,k]
    Nm[1, k] <- B[1, k]
    for (i in 1:(s-1)) {
      B[i+1, k] <- Bstar[i+1, k] * (phi[i,k]-1)/log(phi[i,k])
      Nm[i+1, k] <- Np[i, k] * phi[i, k] + B[i+1, k]
      t_aucn[i,k]<-0.5 * time_AUC[i+1]* (Np[i,k] + Nm[i+1,k])
      t_aucnb[i,k] <- 0.5*(time_AUC[i+1]*Np[i,k]*(phi[i,k]+1))+time_AUC[i+1] *B[i+1, k] *((phi[i,k]/(phi[i,k]-1))-(1/(log(phi[i,k]))))
      }
    for (i in 1:s) {
      Np[i,k] <- Nm[i,k] - n[i,k] +R[i,k]
    }
    AUCn[k] <- 0.5*time_AUC[1]*Nm[1,k]+0.5*time_AUC[s]*Np[s,k]+sum(t_aucn[1:(s-1),k])
    t_r1[k] <- time_AUC[1] *B[1,k]*0.5
    t_r2[k] <- 0.5*time_AUC[s]*Np[s,k]

    AUCc[k] <- t_r1[k] + t_r2[k] + sum(t_aucnb[1:(s-1), k])
    RT[k] <- AUCc[k]/ Nsuper[k]      # equation 1 in Manske & Schwarz 2000
    ESC[k] <- round(AUCn[k]/RT[k])   # equation 2 in Manske & Schwarz 2000
  }
  
####################################################################################################################### #
###########              Partition Bstar estimates by proportional data: race, clips, sex and age             ###########
####################################################################################################################### #
  #derived
  Sigma_pclip<-inverse(inv_Sigma_pclip)
  #Sigma_p_female<-inverse(inv_Sigma_p_female)
  for(i in 1:s){
      ptule[i] <- ilogit(logit_ptule[i])  
  }
  for(k in 1:num_strata){
  sigma_pclip_process[k] <- sqrt(Sigma_pclip[k,k])
  #sigma_p_female_process[k] <- sqrt(Sigma_p_female[k,k])
    for (j in 1:num_strata){
      rho_pclip[k,j] <- (Sigma_pclip[k,j]/(sigma_pclip_process[k]*sigma_pclip_process[j]))
      #rho_p_female[k,j] <- (Sigma_p_female[k,j]/(sigma_p_female_process[k]*sigma_p_female_process[j]))
    }
    for(i in 1:s){
      pclip[i,k] <-ilogit(logit_pclip[i,k])
      #pfemale[i,k] <-ilogit(logit_p_female[i,k])
      for(a in 1:num_ages){
        delta_age[i, k, a] <- exp(log_delta_age[i,k,a])
        p.age[i, k, a] <- delta_age[i,k,a]/sum(delta_age[i,k,1:num_ages])
      }
    }
  }
  #priors
  inv_Sigma_pclip ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  #inv_Sigma_p_female ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  sigma_ptule ~ dt(hyper_value_sigma_ptule_mu,hyper_value_sigma_ptule_sd^-2,1) T(0,)
  sigma_p.age  ~ dt(hyper_value_sigma_p.age_mu,hyper_value_sigma_p.age_sd^-2,1) T(0,)
  for(i in 2:s){
    logit_ptule[i] ~ dnorm(logit_ptule[i-1],sigma_ptule^-2) #p.tule are a univariate timeseries (only one dataset regardless of the total number of strata - not enough data to break up by strata)
    logit_pclip[i,1:num_strata] ~ dmnorm(logit_pclip[i-1,1:num_strata],inv_Sigma_pclip) #p.clip are MVN timeseries among strata
    #logit_p_female[i,1:num_strata] ~ dmnorm(logit_p_female[i-1,1:num_strata],inv_Sigma_p_female) #p.female are MVN timeseries among strata
  }
  
  logit_ptule[1] ~ dnorm(hyper_value_logit_ptule_1_mu, hyper_value_logit_ptule_1_sd^-2) #initiate 1st state for ptule
  for(k in 1:num_strata){
    logit_pclip[1,k] ~ dnorm(hyper_value_logit_pclip_1_mu, hyper_value_logit_pclip_1_sd^-2) #initiate 1st state for each strata
    #logit_p_female[1,k] ~ dnorm(0, 5^-2) #initiate 1st state for each strata
    for(a in 1:num_ages){
      log_delta_age[1, k, a] ~ dnorm(hyper_value_log_delta_age_1_mu,hyper_value_log_delta_age_1_sd^-2) #initiate 1st state for each strata & age
      for(i in 2:s){
        log_delta_age[i, k, a] ~ dnorm(log_delta_age[i-1, k, a],sigma_p.age^-2) #p.ages are univariate timeseries among ages and strata
      }
    }
  }
  #Likelihoods
  for(k in 1:num_strata){
    for(i in 1:s){
      ex_clip[i, k] ~ dbin(pclip[i, k],ex_adfin[i, k])    # proportion adipose clipped
      #ex_female[i, k] ~ dbin(pfemale[i, k],ex_sex[i, k]) # proportion female
      age_dat[i, k, 1:num_ages] ~ dmulti(p.age[i, k, 1:num_ages], age_tot[i, k]) # proportion by age-class
    }
  }
  for(i in 1:s){
      ex_tule[i] ~ dbin(ptule[i], ex_race[i])             # proportion tule (outside k-loop b/c only one grouping of race data)
  }
}
