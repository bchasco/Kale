library(dplyr)
library(tidyr)

d <- read.csv("../data/simpleData2.csv") %>%
  mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
         r_wk = lubridate::week(lubridate::mdy(RecapDate)),
         t_yr = lubridate::year(lubridate::mdy(TagDate)),
         r_yr = lubridate::year(lubridate::mdy(RecapDate))) %>%
  filter(t_yr == 2023) %>%
  filter(t_wk > 10) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  mutate(y_i = as.integer(factor(Return_Yr))) %>%
  filter(is.na(r_wk) | r_k>0) %>%
  mutate(tag = ifelse(Tag1=="",FALSE,TRUE)) %>%
  group_by(y_i,t_k,r_k,t_l, r_l, tag) %>%
  summarise(n = n())

library(RTMB)

# Load Data
data <- list(
  y_i = d$y_i, #year effect
  t_k = d$t_k,   # First detection week
  r_k = d$r_k,   # Recapture week (NA if not recaptured)
  t_l = d$t_l,   # First detection location
  r_l = d$r_l,   # Recapture location (NA if not recaptured)
  tag = as.integer(d$tag), # Whether tagged (1 = tagged, 0 = not)
  n = d$n,
  n_y = length(unique(d$y_i))
  )

parameters = list(
  surv_par = matrix(-2,data$n_y,15),      # Logit persistence probability
  detection_par = matrix(0,data$n_y,5),        # Logit detection probability
  taggingRate_par = matrix(0,data$n_y,5),      # Logit tagging probability
  taggingRate_re = array(0,c(data$n_y,5,max(data$t_k))),      # Logit tagging probability
  par_PopTotal = rep(0,data$n_y),   # Log expected carcass births
  B_time_sig = 0,
  B_loc_sig = 0,
  taggingRate_sig = 0,
  B_time = matrix(0, data$n_y, max(data$t_k)),  #Carcasses by week
  B_loc = array(0,c(data$n_y, 5,max(data$t_k)))  #Distribution of carcasses by location
)

# Define RTMB Model
rtmb_model <- function(parms){

    RTMB::getAll(data, parms)

    nll <- 0
    nll_detect <- array(1,c(n_y,6,max(t_k)))
    nll_tag <- array(0,c(n_y,6,max(t_k)))
    nll_CJS <- rep(0,length(t_k))

    #Survival
    surv <- array(0,c(n_y,6,6))
    surv2 <- array(0,c(n_y,6,6))
    detection <- array(0,c(n_y,6,6))
    #Survival
    for(y in unique(y_i)){
      ii <- 1
      for(i in 1:5){
        for(j in (i+1):6){
          surv[y,i,j] <- exp(surv_par[y,ii])
          ii <- ii + 1
        }
        surv[y,i,i] <- -sum(surv[y,i,])
      }
      surv[y,,] <- as.matrix(Matrix::expm(surv[y,,]))

      surv2[y,,] <- t(surv[y,,])
    }

    # Compute probability matrix
    p <- array(0,c(n_y,6,6))
    for(y in unique(y_i)){
      diag(p[y,,]) <- c(plogis(detection_par[y,]),0)
      p[y,,6] <- 1 - rowSums(p[y,,1:5])
      p[y,6,6] <- 1
    }


    # Derived variables
    N <- array(0, c(n_y,6,max(t_k)))      # Total available carcasses over time
    TotalCarcasses_t <- array(0,c(n_y,6,max(t_k)))      # Observed Total carcasses over time
    TaggedCarcasses_t <- array(0,c(n_y,6,max(t_k)))      # Observed Tagged carcasses over time
    E_TotalCarcasses_t <- array(0,c(n_y,6,max(t_k)))      # Expected Total carcasses over time
    E_TaggedCarcasses_t <- array(0,c(n_y,6,max(t_k)))      # Expected Tagged carcasses over time

    #Conditional probability of loc for each carcass
    B_loc_dist <- array(0,c(n_y,6,max(t_k)))
    for(y in unique(y_i)){
      B_loc_dist[y,,1] <- exp(c(B_loc[y,,1],0))/sum(exp(c(B_loc[y,,1],0)))
      # Initialize spatial distribution for the first time step
      N[y,,1] <- exp(B_time[y,1] + par_PopTotal[y]) * B_loc_dist[y,,1]
      # Loop through time steps to compute derived variables
      for (t in 2:max(t_k)) {

        #Births by time and location
        B_loc_dist[y,,t] <- exp(c(B_loc[y,,t],0))/sum(exp(c(B_loc[y,,t],0)))
        # Total carcasses at time t (new + survivors)
        N[y,,t] = t(surv[y,,]) %*% ((exp(B_time[y,t] + par_PopTotal[y]) * B_loc_dist[y,,t]) +  N[y,,t-1])
      }
    }



    B <- exp(B_time + par_PopTotal)
    nll_B <- rep(0,max(t_k))

    for(y in unique(y_i)){
      for (t in unique(t_k)) {
        # Detection Process (Binomial likelihood)
        for(loc in 1:5){
          TaggedCarcasses_t[y,loc,t] <- sum(na.omit(n[y_i == y & t_k == t & t_l == loc & tag == 1]))  # Tagged carcasses at time t
          TotalCarcasses_t[y,loc,t] <- sum(na.omit(n[y_i == y & t_k == t & t_l == loc]))  # Total detected carcasses at time t
        }

        #Carcass detection probability likelihood
        E_TotalCarcasses_t[y,1:5,t] <- N[y,1:5,t]*plogis(detection_par[y,])
        nll_detect[y,1:5,t] <- RTMB::dpois(TotalCarcasses_t[y,1:5,t]+1,
                                         E_TotalCarcasses_t[y,1:5,t],
                                         log = TRUE)

        #Carcass tagging rate likelihood
        for(loc in 1:5){
          E_TaggedCarcasses_t[y,loc,t] <- TotalCarcasses_t[y,loc,t] * plogis(taggingRate_par[y,loc] + taggingRate_re[y,loc,t])
          if(TotalCarcasses_t[y,loc,t]>0)
            nll_tag[y,loc,t] <- RTMB::dpois(TaggedCarcasses_t[y,loc,t],
                                          (TotalCarcasses_t[y,loc,t] + 0.1) * plogis(taggingRate_par[y,loc] + taggingRate_re[y,loc,t]),
                                          log = TRUE)
        }
      }
    }



    for(i in 1:length(t_l)){
      if(tag[i]==TRUE){
        m <- matrix(0,6,6)
        diag(m) <- 1
        last_k <- ifelse(is.na(r_k[i]),max(na.omit(r_k)),r_k[i])

        if(is.na(r_k[i])){#never recapped
          for(k in (t_k[i]+1):last_k){
            m <- m %*% surv[y_i[i],,] %*% diag(p[y_i[i],,6])
          }
          # m <- m %*% surv %*% diag(p[,6])
        }

        if(!is.na(r_k[i])){#Recaptured
          if(r_k[i]>0){
            # print(r_k[i])
            for(k in (t_k[i]+1):r_k[i]){
              if(k<last_k){
                m <- m %*% surv[y_i[i],,] %*% diag(p[y_i[i],,6]) #non-detection
              }else{
                m <- m %*% surv[y_i[i],,] %*% diag(p[y_i[i],,r_l[i]]) #recap
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

    taggingRate <- taggingRate_re * 0
    B_ts <- B_loc_dist * 0
    B_total <- rep(0,n_y)

    for(y in unique(y_i)){
      taggingRate[y,,] <- plogis(taggingRate_re[y,,] + taggingRate_par[y,])
      B_ts[y,,] <- log(t(B[y,] * t(B_loc_dist[y,,])))
      B_total[y] <- log(sum(B[y,]))
    }

    probability_of_outcome <- nll_CJS

    E <- d %>%
      group_by(t_k, t_l, tag) %>%
      mutate(t_k_total = sum(n))  %>%
      mutate(obs = n)
    E$p <-  probability_of_outcome
    pred2 <- log(E$p * E$t_k_total)


    REPORT(N)
    REPORT(B)
    REPORT(B_ts)
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
sdr <- sdreport(obj)
sd.est <- as.list(sdr,"Estimate", report = TRUE)
sd.sd <- as.list(sdr,"Std. Error", report = TRUE)

