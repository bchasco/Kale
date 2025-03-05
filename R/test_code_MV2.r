library(dplyr)
library(tidyr)

d <- read.csv("../data/simpleData2.csv") %>%
  mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
         r_wk = lubridate::week(lubridate::mdy(RecapDate)),
         t_yr = lubridate::year(lubridate::mdy(TagDate)),
         r_yr = lubridate::year(lubridate::mdy(RecapDate))) %>%
  filter(t_yr == 2022) %>%
  filter(t_wk > 10) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  filter(is.na(r_wk) | r_k>0) %>%
  mutate(tag = ifelse(Tag1=="",FALSE,TRUE)) %>%
  group_by(t_k,r_k,t_l, r_l, tag) %>%
  summarise(n = n())

library(RTMB)

# Load Data
data <- list(
  t_k = d$t_k,   # First detection week
  r_k = d$r_k,   # Recapture week (NA if not recaptured)
  t_l = d$t_l,   # First detection location
  r_l = d$r_l,   # Recapture location (NA if not recaptured)
  tag = as.integer(d$tag), # Whether tagged (1 = tagged, 0 = not)
  n = d$n
  )

parameters = list(
  surv_par = rep(-2,15),      # Logit persistence probability
  detection_par = rep(0,5),        # Logit detection probability
  taggingRate_par = rep(0,5),      # Logit tagging probability
  taggingRate_re = matrix(0,5,max(data$t_k)),      # Logit tagging probability
  par_PopTotal = 1,   # Log expected carcass births
  B_time_sig = 0,
  B_loc_sig = 0,
  taggingRate_sig = 0,
  B_time = rep(0, max(data$t_k)),  #Carcasses by week
  B_loc = matrix(0, 5,max(data$t_k))  #Distribution of carcasses by location
)

# Define RTMB Model
rtmb_model <- function(parms){

    RTMB::getAll(data, parms)

    nll <- 0
    nll_detect <- matrix(1,6,max(t_k))
    nll_tag <- matrix(0,6,max(t_k))
    nll_CJS <- rep(0,length(t_k))

    #Survival
    surv <- matrix(0,6,6)
    detection <- matrix(0,6,6)
    #Survival
    ii <- 1
    for(i in 1:5){
      for(j in (i+1):6){
        surv[i,j] <- exp(surv_par[ii])
        ii <- ii + 1
      }
      surv[i,i] <- -sum(surv[i,])
    }
    surv <- as.matrix(Matrix::expm(surv))

    surv2 <- t(surv)

    # Compute probability matrix
    p <- matrix(0,6,6)
    diag(p) <- c(plogis(detection_par),0)
    p[,6] <- 1 - rowSums(p[,1:5])
    p[6,6] <- 1


    # Derived variables
    N <- matrix(0,6,max(t_k))      # Total available carcasses over time
    TotalCarcasses_t <- matrix(0,6,max(t_k))      # Observed Total carcasses over time
    TaggedCarcasses_t <- matrix(0,6,max(t_k))      # Observed Tagged carcasses over time
    E_TotalCarcasses_t <- matrix(0,6,max(t_k))      # Expected Total carcasses over time
    E_TaggedCarcasses_t <- matrix(0,6,max(t_k))      # Expected Tagged carcasses over time

    #Conditional probability of loc for each carcass
    B_loc_dist <- matrix(0,6,max(t_k))
    B_loc_dist[,1] <- exp(c(B_loc[,1],0))/sum(exp(c(B_loc[,1],0)))

    # Initialize spatial distribution for the first time step
    N[,1] <- exp(B_time[1] + par_PopTotal) * B_loc_dist[,1]

    # Loop through time steps to compute derived variables
    for (t in 2:max(t_k)) {

      #Births by time and location
      B_loc_dist[,t] <- exp(c(B_loc[,t],0))/sum(exp(c(B_loc[,t],0)))
      # Total carcasses at time t (new + survivors)
      N[,t] = t(surv) %*% ((exp(B_time[t] + par_PopTotal) * B_loc_dist[,t]) +  N[,t-1])
    }

    B <- exp(B_time + par_PopTotal)
    nll_B <- rep(0,max(t_k))

    for (t in unique(t_k)) {
      # Detection Process (Binomial likelihood)
      for(loc in 1:5){
        TaggedCarcasses_t[loc,t] <- sum(na.omit(n[t_k == t & t_l == loc & tag == 1]))  # Tagged carcasses at time t
        TotalCarcasses_t[loc,t] <- sum(na.omit(n[t_k == t & t_l == loc]))  # Total detected carcasses at time t
      }

      #Carcass detection probability likelihood
      E_TotalCarcasses_t[1:5,t] <- N[1:5,t]*plogis(detection_par)
      nll_detect[1:5,t] <- RTMB::dpois(TotalCarcasses_t[1:5,t],
                                       E_TotalCarcasses_t[1:5,t],
                                       log = TRUE)

      #Carcass tagging rate likelihood
      for(loc in 1:5){
        E_TaggedCarcasses_t[loc,t] <- TotalCarcasses_t[loc,t] * plogis(taggingRate_par[loc] + taggingRate_re[loc,t])
        if(TotalCarcasses_t[loc,t]>0)
        nll_tag[loc,t] <- RTMB::dpois(TaggedCarcasses_t[loc,t],
                                       TotalCarcasses_t[loc,t] * plogis(taggingRate_par[loc] + taggingRate_re[loc,t]),
                                       log = TRUE)
      }
    }



    for(i in 1:length(t_l)){
      if(tag[i]==TRUE){
        m <- matrix(0,6,6)
        diag(m) <- 1
        last_k <- ifelse(is.na(r_k[i]),max(na.omit(r_k)),r_k[i])

        if(is.na(r_k[i])){#never recapped
          for(k in (t_k[i]+1):last_k){
            m <- m %*% surv %*% diag(p[,6])
          }
          # m <- m %*% surv %*% diag(p[,6])
        }

        if(!is.na(r_k[i])){#Recaptured
          if(r_k[i]>0){
            # print(r_k[i])
            for(k in (t_k[i]+1):r_k[i]){
              if(k<last_k){
                m <- m %*% surv %*% diag(p[,6]) #non-detection
              }else{
                m <- m %*% surv %*% diag(p[,r_l[i]]) #recap
              }
            }
            last_p <- matrix(0,6,6)
            last_p[6,6] <- 1
            # m <- m %*% surv %*% last_p  #Death
          }
        }
        delta <- rep(0,6)
        delta[t_l[i]] <- 1
        nll_CJS[i] <- t(delta) %*% m %*% rep(1,6)
      }
    }

    nll_CJS_total <- sum(log(nll_CJS[tag==TRUE])*n[tag==TRUE])

    nll_B_time <- RTMB::dnorm(B_time, 0, exp(B_time_sig), log = TRUE)
    nll_B_loc <- RTMB::dnorm(B_loc, 0, exp(B_loc_sig), log = TRUE)
    nll_taggingRate_re <- RTMB::dnorm(taggingRate_re, 0, exp(taggingRate_sig), log = TRUE)

    taggingRate <- plogis(taggingRate_par + taggingRate_re)

    probability_of_outcome <- nll_CJS

    E <- d %>%
      group_by(t_k, t_l, tag) %>%
      mutate(t_k_total = sum(n))  %>%
      mutate(obs = n)
    E$p <-  probability_of_outcome
    pred2 <- log(E$p * E$t_k_total)

    B_ts <- log(t(B * t(B_loc_dist)))

    B_total <- log(sum(B))

    REPORT(N)
    REPORT(B)
    REPORT(B_loc_dist)
    REPORT(surv)
    REPORT(surv2)
    REPORT(p)
    REPORT(probability_of_outcome)
    REPORT(taggingRate)
    REPORT(TaggedCarcasses_t)
    REPORT(TotalCarcasses_t)
    REPORT(E_TaggedCarcasses_t)
    REPORT(E_TotalCarcasses_t)
    REPORT(nll_tag)
    REPORT(nll_detect)
    REPORT(nll_CJS)
    RTMB::ADREPORT(pred2)
    RTMB::ADREPORT(B)
    RTMB::ADREPORT(E_TaggedCarcasses_t)
    RTMB::ADREPORT(E_TotalCarcasses_t)
    RTMB::ADREPORT(B_ts)
    RTMB::ADREPORT(B_total)
    RTMB::ADREPORT(taggingRate)

    return(-(sum(nll_B) +
            nll_CJS_total +
            sum(nll_detect) +
            sum(nll_tag) +
            sum(nll_B_time) +
            sum(nll_B_loc) +
            sum(nll_taggingRate_re))
    )
}

obj <- RTMB::MakeADFun(rtmb_model,
                       parameters,
                       random = c("B_time","B_loc","taggingRate_re"),
                       silent = TRUE)

# Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- obj$report()
sd.est <- as.list(sdreport(obj),"Estimate", report = TRUE)
sd.sd <- as.list(sdreport(obj),"Std. Error", report = TRUE)

