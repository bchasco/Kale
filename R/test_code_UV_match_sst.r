library(fs)
library(tibble)
library(knitr)
library(dplyr)
library(tidyr)
library(readr)


combined_data <- read.csv("./data/NF_Lewis_combined-2013_2024-2025-03-19.csv")
cols <- c("Return_Yr","SPECIES","Run","ScaleAge", "Sex","TagDate","TagReach",
          "Tag1","Tag2","RecapDate","RecapReach")

result <- combined_data %>%
  mutate(oa = ifelse(is.na(ScaleAge)|(ScaleAge==9)|(ScaleAge==88),-1,as.numeric(substr(ScaleAge,1,1)) - as.numeric(substr(ScaleAge,2,2)))) %>%
  mutate(Sex = ifelse(is.na(Sex)|(Sex=="Unknown"),-1,Sex)) %>%
  mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
         r_wk = ifelse(is.na(RecapDate),NA,lubridate::week(lubridate::dmy(RecapDate))),
         y_i = Return_Yr) %>%
  mutate(r_wk = ifelse(as.numeric(as.character(r_wk))<20,r_wk+ 52,r_wk),
         t_wk = ifelse(as.numeric(as.character(t_wk))<20,t_wk+ 52,t_wk)) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1) %>%
  group_by(t_wk) %>%

# Create lookup tables separately for TagReach and RecapReach
lk_up <- data.frame(
  TagReach = sort(unique(result$TagReach)),
  t_l = seq_along(unique(result$TagReach))
)

lk_up_r <- data.frame(
  TagReach = sort(unique(result$TagReach)),
  r_l = seq_along(unique(result$TagReach))
)
#
# Perform left joins to assign t_l and r_l separately
year <- 2013
d <- result %>%
  filter(t_k<=22
         ,Return_Yr %in% c(!!year)
         ) %>%
  left_join(lk_up, by = c("TagReach" = "TagReach")) %>%  # Adds t_l based on TagReach
  left_join(lk_up_r, by = c("RecapReach" = "TagReach")) %>%   # Adds r_l based on RecapReach
  mutate(check = ifelse(is.na(RecapReach),1,ifelse(RecapReach<TagReach,2,1))) %>%
  filter(check==1) %>%
  mutate(t_l = 1,
         r_l = ifelse(is.na(r_l),NA,1)) %>%
  mutate(tag = ifelse(is.na(Tag1),FALSE,TRUE)) %>%
  mutate(t_k = t_k - min(t_k) + 1,
         r_k = r_k - min(t_k) + 1) %>%
  group_by(y_i,t_k,r_k,TagReach,RecapReach,t_l, r_l, tag) %>%
  summarise(n = n())

library(RTMB)

max_l <- max(na.omit(d$r_l)) + 1

# Load Data
data <- list(
  y_i = as.integer(as.factor(d$y_i)), #year effect
  t_k = d$t_k,   # First detection week
  r_k = d$r_k,   # Recapture week (NA if not recaptured)
  t_l = d$t_l,   # First detection location
  r_l = d$r_l,   # Recapture location (NA if not recaptured)
  tag = as.integer(d$tag), # Whether tagged (1 = tagged, 0 = not)
  n = d$n,
  n_y = length(unique(d$y_i))
  )

parameters = list(
  surv_par = 0,      # Logit persistence probability
  detect_par = 0,
  taggingRate_par = 0,      # Logit tagging probability
  par_PopTotal = 6,   # Log expected carcass births
  B_time_sig = 0,
  B_time = rep(0, max(d$t_k))
)
#
# # Define RTMB Model
rtmb_model <- function(parms){

    RTMB::getAll(data, parms)

    nll <- 0

    max_l <- max(na.omit(r_l)) + 1

    print(max_l)

    nll_detect <- rep(1,max(t_k))
    nll_tag <- rep(0,max(t_k))
    nll_CJS <- rep(0,length(t_k))

    #Survival
    surv <- matrix(0,2,2)
    f_surv <- matrix(0,2,2)
    #Survival
    ii <- 1
    surv[1,2] <- exp(surv_par)
    surv[1,1] <- -surv[1,2]
    f_surv <- surv
    surv <- as.matrix(Matrix::expm(surv))

    # Compute probability matrix
    p <- matrix(0,2,2)
    p[1,2] <- exp(detect_par)
    p[1,1] <- -p[1,2]
    p <- as.matrix(Matrix::expm(p))

    # Derived variables
    N <- matrix(0, 2, max(t_k))      # Total available carcasses over time

    # Initialize spatial distribution for the first time step
    N[,1] <- t(as.matrix(Matrix::expm(f_surv * 0.5))) %*% c(exp(B_time[1] + par_PopTotal),0)
    # Loop through time steps to compute derived variables
    for(t in 2:max(t_k)){
      N[,t] = t(as.matrix(Matrix::expm(f_surv * 0.5))) %*% c(exp(B_time[t] + par_PopTotal),0) +
        t(surv) %*% N[,t-1]
    }


    TotalCarcasses_t <- rep(0,max(t_k)) #Total
    CarcassesW_tags_t <- rep(0,max(t_k)) #tags
    E_TotalCarcasses_t <- matrix(0,2,max(t_k)) #predicted
    E_TaggedCarcasses_t <- matrix(0,max(t_k)) #predicted
    for (t in 1:max(t_k)) {
      # Detection Process (Binomial likelihood)
      CarcassesW_tags_t[t] <- sum(na.omit(n[t_k == t & t_l == 1 & tag == 1]))  # Tagged carcasses at time t
      TotalCarcasses_t[t] <- sum(na.omit(n[t_k == t & t_l == 1]))  # Total detected carcasses at time t

      #Carcass detection probability likelihood
      E_TotalCarcasses_t[,t] <- t(p) %*% as.vector(N[,t])
      nll_detect[t] <- RTMB::dpois(TotalCarcasses_t[t]+1,
                                   t(c(1,0)) %*% E_TotalCarcasses_t[,t],
                                       log = TRUE)
      #Carcass tagging rate likelihood
      E_TaggedCarcasses_t[t] <- E_TotalCarcasses_t[1,t] * plogis(taggingRate_par)

      nll_tag[t] <- RTMB::dpois(CarcassesW_tags_t[t] + 1,
                                      E_TaggedCarcasses_t[t],
                                    log = TRUE)
      # }
    }
    REPORT(TotalCarcasses_t)
    REPORT(E_TotalCarcasses_t)
    REPORT(E_TaggedCarcasses_t)
    REPORT(nll_tag)
    REPORT(nll_detect)
    REPORT(N)

    for(i in 1:length(t_l)){
      m <- matrix(0,max_l,max_l)
      diag(m) <- 1
      last_k <- ifelse(is.na(r_k[i]),max(na.omit(r_k)),r_k[i])

      if(is.na(r_k[i])){#never recapped
        if((t_k[i]+1)<last_k){
          for(k in (t_k[i]+1):last_k){
            m <- m %*% surv %*% diag(p[,max_l])
          }
        }
      }


      if(!is.na(r_k[i])){#Recaptured
        if(r_k[i]>0){
          for(k in (t_k[i]+1):r_k[i]){
            if(k<last_k){
              m <- m %*% surv %*% diag(p[,max_l]) #non-detection
            }else{
              m <- m %*% surv %*% diag(p[,r_l[i]]) #recap
            }
          }
          last_p <- matrix(0,max_l,max_l)
          last_p[max_l,max_l] <- 1
          # m <- m %*% surv %*% last_p  #Death
        }
      }
      delta <- rep(0,max_l)
      delta[t_l[i]] <- 1
      nll_CJS[i] <- t(delta) %*% m %*% rep(1,max_l)
    }

    nll_CJS_total <- sum(log(nll_CJS[tag==TRUE])*n[tag==TRUE])
    nll_B_time <- RTMB::dnorm(B_time, -0.5*exp(B_time_sig)^2, exp(B_time_sig), log = TRUE)

    B <- exp(B_time + par_PopTotal)
    Btotal <- sum(B)
    REPORT(B)
    REPORT(surv)
    REPORT(p)
    REPORT(Btotal)
    REPORT(nll_CJS_total)
    REPORT(taggingRate_par)
    ADREPORT(Btotal)
    ADREPORT(surv)
    ADREPORT(E_TotalCarcasses_t)
    ADREPORT(B)
    return(-(nll_CJS_total +
            sum(nll_detect) +
            sum(nll_tag) +
            sum(nll_B_time))
    )
    # return(0)
}
#
obj <- RTMB::MakeADFun(rtmb_model,
                       parameters,
                       random = c("B_time"),
                       silent = FALSE)

# Optimize Model
opt <- TMBhelper::fit_tmb(obj)
rep <- obj$report()
sdr <- sdreport(obj)
as.list(sdr,"Estimate",report=TRUE)$B

# print(round(rep$E_TotalCarcasses_t))
# print(round(rep$N))
# print(year)
# print(sum(rep$B))
# print(sdr$sd)
#
