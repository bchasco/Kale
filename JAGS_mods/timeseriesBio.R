####################################################################################################################### #
###########              Partition Bstar estimates by proportional data: race, clips, sex and age             ###########
####################################################################################################################### #
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
  # for(i in 1:num_ages){
  #   alpha[i]<-1
  # }
}
model{
  #derived
  Sigma_p_tule<-inverse(inv_Sigma_p_tule)
  Sigma_p_clip<-inverse(inv_Sigma_p_clip)
  Sigma_p_female<-inverse(inv_Sigma_p_female)
  for(k in 1:num_strata){
  sigma_p_tule_process[k] <- sqrt(Sigma_p_tule[k,k])
  sigma_p_clip_process[k] <- sqrt(Sigma_p_clip[k,k])
  sigma_p_female_process[k] <- sqrt(Sigma_p_female[k,k])
    for (j in 1:num_strata){
      rho_p_tule[k,j] <- (Sigma_p_tule[k,j]/(sigma_p_tule_process[k]*sigma_p_tule_process[j]))
      rho_p_clip[k,j] <- (Sigma_p_clip[k,j]/(sigma_p_clip_process[k]*sigma_p_clip_process[j]))
      rho_p_female[k,j] <- (Sigma_p_female[k,j]/(sigma_p_female_process[k]*sigma_p_female_process[j]))
    }
    for(i in 1:s){
      ptule[i,k] <- ilogit(logit_p_tule[i,k])
      pclip[i,k] <-ilogit(logit_p_clip[i,k])                               
      pfemale[i,k] <-ilogit(logit_p_female[i,k])     
      for(a in 1:num_ages){
        delta_age[i, k, a] <- exp(log_delta_age[i,k,a])
        p.age[i, k, a] <- delta_age[i,k,a]/sum(delta_age[i,k,1:num_ages])
      }
    }
  }
  #priors
  inv_Sigma_p_tule ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  inv_Sigma_p_clip ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  inv_Sigma_p_female ~ dwish(Rmat[1:num_strata,1:num_strata],num_strata + 1)
  sigma_p_age ~ dt(0,5^-2,1) T(0,)
  for(i in 2:s){
    logit_p_tule[i,1:num_strata] ~ dmnorm(logit_p_tule[i-1,1:num_strata],inv_Sigma_p_tule) #p.tule are MVN timeseries among strata
    logit_p_clip[i,1:num_strata] ~ dmnorm(logit_p_clip[i-1,1:num_strata],inv_Sigma_p_clip) #p.clip are MVN timeseries among strata
    logit_p_female[i,1:num_strata] ~ dmnorm(logit_p_female[i-1,1:num_strata],inv_Sigma_p_female) #p.female are MVN timeseries among strata
  }
  for(k in 1:num_strata){
    logit_p_tule[1,k] ~ dnorm(0, 5^-2) #initiate 1st state for each strata
    logit_p_clip[1,k] ~ dnorm(0, 5^-2) #initiate 1st state for each strata
    logit_p_female[1,k] ~ dnorm(0, 5^-2) #initiate 1st state for each strata
    for(a in 1:num_ages){
      log_delta_age[1, k, a] ~ dnorm(0,3^-2) #initiate 1st state for each strata & age
      for(i in 2:s){
        log_delta_age[i, k, a] ~ dnorm(log_delta_age[i-1, k, a],sigma_p_age^-2) #p.ages are univariate timeseries among ages and strata
      }
    }
  }
  #Likelihoods
  for(k in 1:num_strata){
    for(i in 1:s){
      ex_tule[i, k] ~ dbin(ptule[i, k], ex_race[i, k])               # proportion tule
      ex_clip[i, k] ~ dbin(pclip[i, k],ex_adfin[i, k])               # proportion adipose clipped
      ex_female[i, k] ~ dbin(pfemale[i, k],ex_sex[i, k])           # proportion female
      age_dat[i, k, 1:num_ages] ~ dmulti(p.age[i, k, 1:num_ages], age.tot[i, k]) # proportion by age-class
    }
  }
}