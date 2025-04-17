model{
  for (i in 1:(s-1)){
    R[i] ~dbin(v[i],n[i])
    r[i] ~dbin(lambda[i],R[i])
    psi[i+1] <- psi[i]*(1-p)*phi + pent[i+1] *(phi-1)/log(phi) 
    lambda[i] <- phi*(p+(1-p)*lambda[i+1])
  }
  for (i in 2:(s-1)){
    m[i] ~dbin(tau[i],T[i])
  }
  for (i in 1:s){
    v[i] ~dbeta(a,b)
    delta[i]~dgamma(alpha[i],1)
    pent[i]<-delta[i]/sum(delta[1:s])
    temp[i] <- psi[i]*p
    multP[i] <- temp[i]/sum(temp[1:s])
    tau[i] <-p/(p+(1-p)*lambda[i]) 
  }
  psi[1]<-pent[1]
  phi~dbeta(a,b)
  p~dbeta(a,b)
  lambda[s] <- 0
  Ntot ~ dunif(LL,UL) 
  Nsuper <- round(Ntot);
  u[1:s] ~ dmulti(multP[1:s],uTot)
  psiPtot <- sum(temp[1:s])
  uTot ~ dbin(psiPtot,Nsuper)
}
