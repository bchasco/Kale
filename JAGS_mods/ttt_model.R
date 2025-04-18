model{
  for (i in 1:(s-1)){
    phi[i] ~dbeta(a,b)
    R[i] ~dbin(v[i],n[i])
    r[i] ~dbin(lambda[i],R[i])
    psi[i+1] <- psi[i]*(1-p[i])*phi[i] + pent[i+1] *(phi[i]-1)/log(phi[i]) 
    lambda[i] <- phi[i]*(p[i+1]+(1-p[i+1])*lambda[i+1])
  }
  for (i in 2:(s-1)){
    p[i] ~dbeta(a,b)
    m[i] ~dbin(tau[i],T[i])
  }
  for (i in 1:s){
    v[i] ~dbeta(a,b)
    delta[i]~dgamma(alpha[i],1)
    pent[i]<-delta[i]/sum(delta[1:s])
    temp[i] <- psi[i]*p[i]
    multP[i] <- temp[i]/sum(temp[1:s])
    tau[i] <-p[i]/(p[i]+(1-p[i])*lambda[i]) 
  } 
  psi[1]<-pent[1]
  p[1]<-p[2]
  p[s]<-p[s-1]
  lambda[s] <- 0
  Ntot ~ dunif(LL,UL) 
  Nsuper <- round(Ntot);
  u[1:s] ~ dmulti(multP[1:s],uTot)
  psiPtot <- sum(temp[1:s])
  uTot ~ dbin(psiPtot,Nsuper)
}
  