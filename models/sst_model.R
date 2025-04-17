model{

  for (i in 1:s){
    v[i] ~dbeta(a,b) #prior for marking rate
    delta[i]~dgamma(alpha[i],1) #NUMBER of carcasses entering
    pent[i]<-delta[i]/sum(delta[1:s]) #scaled probability of entry
    temp[i] <- psi[i]*p #joint survival after entry and then detected
    multP[i] <- temp[i]/sum(temp[1:s])
    tau[i] <-p/(p+(1-p)*lambda[i])  #Recapture probabiilty
  }

  for (i in 1:(s-1)){
    #i is period
    psi[i+1] <- psi[i]*(1-p)*phi + pent[i+1] *(phi-1)/log(phi)  #survival and still availabe for capture
    lambda[i] <- phi*(p+(1-p)*lambda[i+1]) #survival * detection probability * lambda (chi)
    R[i] ~dbin(v[i],n[i]) #marking rate; v = total marks(R[i])/total(mark + unmarked)
    r[i] ~dbin(lambda[i],R[i]) #recapture rate, lambda = total recaps (r)/ total marks (R[i])
  }
  for (i in 2:(s-1)){
    m[i] ~dbin(tau[i],T[i]) #tau = number of previously marked fish (m[i]) / marks to recapture (T = m + z)
  }


  psi[1]<-pent[1]
  pent_tmp[1]<-pent[1]
  psi_tmp[1]<-pent[1]
  phi~dbeta(a,b)
  p~dbeta(a,b)
  lambda[s] <- 0
  Ntot ~ dunif(LL,UL)
  Nsuper <- round(Ntot);
  u[1:s] ~ dmulti(multP[1:s],uTot)
  est_r <- lambda * R
  est_R <- v * n
  est_u <- multP * uTot
  psiPtot <- sum(temp[1:s])
  est_uTot <- psiPtot*Nsuper
  uTot ~ dpois(est_uTot)
}
