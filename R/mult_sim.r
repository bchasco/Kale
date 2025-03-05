library(RTMB)

nT <- 3 #number of time steps
n <- 100000 #number of individuals

ch <- matrix(0,n,nT) #capture history
true_state <- matrix(0,n,nT) #Simulated true state of the fish
ns <- 2 #number states

#create the survival matrix
surv <- matrix(0,ns,ns)
surv[1,1] <- -0.3
surv[1,2] <- -surv[1,1]
surv <- as.matrix(Matrix::expm(surv))

#Create the detection matrix
p <- matrix(0,ns,ns)
p[1,1] <- -2
p[1,2] <- -p[1,1]
p <- as.matrix(Matrix::expm(p))

#Simulate data
for(i in 1:n){
  s <- matrix(0,ns,nT)
  s[1,1] <- 1 #initial condition
  for(tt in 2:nT){
    s[,tt] <- rmultinom(1,1,(t(surv) %*% s[,tt-1]))
    true_state[i,tt] <- which(s[,tt]==1,arr.ind=TRUE)
    tmp_ch <- rmultinom(1,1,p[true_state[i,tt],])
    ch[i,tt] <- which(tmp_ch==1,arr.ind=TRUE)[1,1]
  }
}

#Parameter and data lists for RTMB
parameters <- list(par = c(0,0))
data <- list(ch = ch)

#Create the model function to be passed to AD libraries.
f <- function(parms){

  RTMB::getAll(data, parms)

  ns <- 2
  surv <- matrix(0,ns,ns)
  surv[1,1] <- -exp(par[1])
  surv[1,2] <- -surv[1,1]
  surv <- as.matrix(Matrix::expm(surv))

  p <- matrix(0,ns,ns)
  p[1,1] <- -exp(par[2])
  p[1,2] <- -p[1,1]
  p <- as.matrix(Matrix::expm(p))

  nll_i <- rep(0, nrow(ch))


  # Compute the negative log-likelihood
  #This is a very, very compact way of writing the generalized likelihood
  for(i in 1:nrow(ch)){
    m <- matrix(0,ns,ns)
    diag(m) <- 1
    for(tt in 2:ncol(ch)){
        m <- m %*% surv %*% diag(p[,ch[i,tt]])
    }
    delta <- rep(0,ns)
    delta[1] <- 1
    ones <- rep(1,ns)
    nll_i[i] <- t(delta) %*% m %*% ones
  }

  RTMB::REPORT(surv)
  RTMB::REPORT(p)
  return(-sum(log(nll_i)))

}

obj <- RTMB::MakeADFun(f,
                       parameters,
                       silent = FALSE)

st <- Sys.time()
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- obj$report()
#Predicted
print("Predicted survival")
print(round(rep$surv,2))
print("Simulated values")
print(round(surv,2))
print("Predicted survival")
print(round(rep$p,2))
print("Simulated values")
print(round(p,2))

print(paste("Total time to estimate ", n, "observations"))
print(Sys.time()-st)
