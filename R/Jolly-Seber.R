library(expm)  # Matrix exponentiation
library(tidyr)
library(dplyr)

rm(list=ls())

#Read in data
d <- read.csv("data/simpleData2.csv") %>%
  mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
         r_wk = lubridate::week(lubridate::mdy(RecapDate))) %>%
  filter(t_wk > 10) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  filter(!(Tag1==""))


data <- list(t_l = d$t_l, #tagging location
             r_l = d$r_l, #recapture location
             t_k = d$t_k, #tagging week
             r_k = d$r_k) #recapture week, last week if not recapture

# Initial parameter values
parameters <- list(
  phi_par = rep(0,15),
  p_par = rep(-1,5)
)

f <- function(parms){

  RTMB::getAll(data,
               parms)

  nll <- rep(0,length(t_l))

  #Survival
  phi <- matrix(0,6,6)
  ii <- 1
  for(i in 1:5){
    for(j in i:5){
      phi[i,j] <- exp(phi_par[ii])
      ii <- ii + 1
    }
    phi[,6] <- 1
  }
  for(i in 1:5){
    phi[i,] <- phi[i,]/sum(phi[i,])
  }

  #Detection probability
  p <- matrix(0,6,6)
  for(i in 1:5) p[i,i] <- RTMB::plogis(p_par[i])
  p[,6] <- c(1-RTMB::plogis(p_par),1)

  # Compute probability
  for(i in 1:length(t_l)){
    m <- matrix(0,6,6)
    diag(m) <- 1
    last_k <- ifelse(is.na(r_k[i]),max(na.omit(r_k)),r_k[i])
    if(is.na(r_k[i])){#never recapped
      for(k in (t_k[i]+1):last_k){
        m <- m %*% phi %*% diag(p[,6])
      }
    }else{
      for(k in (t_k[i]+1):r_k[i]){
        if(k<last_k){
          m <- m %*% phi %*% diag(p[,6]) #non-detection
        }else{
          m <- m %*% (phi) %*% diag(p[,r_l[i]]) #recap
        }
      }
      # m <- m %*% diag(phi[,6]) #Death
    }
    delta <- rep(0,6)
    delta[t_l[i]] <- 1
    nll[i] <- t(delta) %*% m %*% rep(1,6)
  }

  RTMB::REPORT(phi)
  RTMB::REPORT(p)
  RTMB::REPORT(nll)

  q_phi <- RTMB::qlogis(phi) #Survival
  q_p <- RTMB::qlogis(p) #detection probability

  q_phi_p <- RTMB::qlogis((phi%*%p)[1:5,1:6])
  indices <- which(upper.tri(q_phi_p, diag=FALSE), arr.ind=TRUE)
  # mat[indices[order(indices[,1]),]]
  RTMB::REPORT(q_phi)
  RTMB::REPORT(q_p)
  RTMB::REPORT(q_phi_p)
  RTMB::ADREPORT(q_phi_p)

  return(-sum(log(nll)))
}

obj <- RTMB::MakeADFun(f,
                       parameters,
                       silent = FALSE)
opt <- nlminb(obj$par,
              obj$fn,
              obj$gr)


rep <- obj$report()
sd <- RTMB::sdreport(obj)

source("R/plot_the_data.r")
