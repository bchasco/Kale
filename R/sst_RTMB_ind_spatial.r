#---------------------------------------------------------------------------------------------------------- -
#                                                   ----
#---------------------------------------------------------------------------------------------------------- -
# install_or_load_pack <- function(packages) {
#   new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new_packages)) {
#     install.packages(new_packages, dependencies = TRUE)
#   }
#   sapply(packages, require, character.only = TRUE)
# }

# Load functions
# sapply(FUN = source, paste(getwd(), "functions", list.files("functions"), sep="/"))
logit<-function(x){log(x/(1-x))}

# Install/Load packages
package_list<-c("here", "gplots", "tidyverse", "dataRetrieval", "RColorBrewer", "openxlsx"
                , "RMark" , "R2jags", "tidybayes", "MCMCvis", "patchwork", "ggthemes", "glue", "coda")

sapply(package_list, require, character.only = TRUE)

# install_or_load_pack(package_list)

# Specify model files of interest
models<-c(
  "sst_model.r"
)

# Load M-R dataset via saved .rds file
location<-c("NFLewis")
species<-c("Chinook")

JS_stats_ALL<-readRDS(paste0("Data/Originals/NF_Lewis/Lewis_",year,"_JS_data.rds"))

# Collapse 3rd dimension of JS_stats into a single grouping using "collapse_third_dimension" function
# Function to collapse third dimension
collapse_third_dimension <- function(data_list) {
  # Initialize an empty list to store long format data frames
  long_dfs <- list()
  # Loop over each element in the list
  for (level_name in names(data_list)) {
    long_df <- data_list[[level_name]] %>%
      mutate(type = level_name) %>%
      pivot_longer(-c(dates, type), names_to = 'variable', values_to = 'value') %>%
      mutate(value = as.numeric(as.character(value)))
    long_dfs[[level_name]] <- long_df
  }
  # Combine all long data frames
  combined_long <- bind_rows(long_dfs)
  # Summarize to collapse the levels
  collapsed_df <- combined_long %>%
    group_by(dates, variable) %>%
    summarize(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = variable, values_from = value) |>
    mutate(dates = as.Date(dates))
  return(collapsed_df)
}
# Collapse 3rd dimension of JS_stats into a single grouping using "collapse_third_dimension" function
(collapsed_df <- collapse_third_dimension(JS_stats_ALL))
# collapsed_df <- (JS_stats_ALL$Male)
# for(i in 2:7)
#   collapsed_df[,i] <- as.numeric(as.character(collapsed_df[,i]))

(year<-year(min(collapsed_df$dates)))
(collapsed_df <- collapse_third_dimension(JS_stats_ALL))

(year<-year(min(collapsed_df$dates)))

# Define "good" JS vectors of data based on rule-sets
good_s_T_mat_1 = apply(as.matrix(collapsed_df$zi + collapsed_df$mi),2,function(x) c(which(x>0),rep(NA,length(which(x==0)))))
good_s_T_mat_2 = apply(apply(good_s_T_mat_1,2,function(x) ifelse(x%in%c(2:(nrow(collapsed_df)-1)),x,NA)),2,function(x) sort(x,na.last=TRUE))
good_s_T_final = apply(good_s_T_mat_2,2,function(x) length(x[!is.na(x)]))

good_s_R_mat_1<-apply(as.matrix(collapsed_df$Ri),2,function(x) c(which(x>0),rep(NA,length(which(x==0)))))
good_s_R_mat_2=apply(apply(good_s_R_mat_1,2,function(x) ifelse(x%in%c(1:(nrow(collapsed_df)-1)),x,NA)),2,function(x) sort(x,na.last=TRUE))
good_s_R_final<-apply(good_s_R_mat_2,2,function(x) length(x[!is.na(x)]))

# d <- read.csv("data/simpleData2.csv") %>%
#   mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
#          r_wk = lubridate::week(lubridate::mdy(RecapDate)),
#          t_yr = lubridate::year(lubridate::mdy(TagDate)),
#          r_yr = lubridate::year(lubridate::mdy(RecapDate))) %>%
#   filter(t_yr == !!year) %>%
#   filter(t_wk > 10) %>%
#   mutate(t_k = as.integer(as.factor(t_wk - min(t_wk) + 1)),
#          r_k = (as.integer(r_wk) - min(t_wk) + 1)) %>%
#   mutate(r_k = r_k - (max(na.omit(r_k) - max(na.omit(t_k))))) %>%
#   mutate(t_l = TagState,
#          r_l = RecapState) %>%
#   filter(is.na(r_wk) | r_k>0) %>%
#   mutate(tag = ifelse(Tag1=="",FALSE,TRUE)) %>%
#   group_by(t_k,r_k,t_l, r_l, tag, t_yr) %>%
#   summarise(n = n())

d <- read.csv("data/NF_Lewis_combined-2013_2024-2025-04-18.csv") %>%
  mutate(t_k = Period_Cap,
         r_k = Period_Recap,
         t_yr = Return_Yr) %>%
  filter(Return_Yr == !!year) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  # filter(is.na(r_wk) | r_k>0) %>%
  mutate(tag = ifelse(is.na(Tag1),FALSE,TRUE)) %>%
  group_by(t_k,r_k,t_l, r_l, tag, t_yr) %>%
  summarise(n = n())



if(is.null(dim(MR_data$R)[2])==TRUE){
  MR_data[lapply(MR_data,class)=="matrix"]<-lapply(MR_data[lapply(MR_data,class)=="matrix"],as.vector)
  MR_data$alpha<-rep(1,MR_data$s)
}

# Create inits for uTot
# uTot<-MR_data$uTot
# inits_func<-function(){
#   inits<-list(
#     "Ntot" = uTot * runif(1,1,5)
#   )
# }

u <- rep(0,max(d$t_k))
u_i <- aggregate(list(u = d$n),by = list(s = d$t_k),sum)
u[u_i$s] <- u_i$u

data <- list(
  s = max(d$t_k),                    # Number of periods
  t_k = d$t_k,               # First detection
  r_k = d$r_k,               # Recapture time (NA if none)
  tag = as.integer(d$tag),   # 1 = tagged, 0 = untagged
  n = d$n,                   # Number of fish in each group
  u = u,          # First detections per period (for multinomial)
  uTot = sum(d$n)           # Total number of detections
  # # same delta setup as before
  # ...
)

print(data$uTot)

## Parameters
parameters <- list(
  f_alpha = log(0.5),
  f_va= qlogis(0.5),
  f_vb= qlogis(0.5),
  logit_phi = 0,        # logit of survival
  logit_p = 0,          # logit of detection
  log_Ntot = 10,        # log total N
  log_delta = rep(-0, data$s),# log of delta (entry intensity)
  logit_v = rep(0, data$s)   # logit marking rate v[i]
)

require(RTMB)

model <- function(parms){

  RTMB::getAll(data,parms)
  ## Transformed parameters
    phi <- plogis(logit_phi)
    p <- plogis(logit_p)
    Ntot <- exp(log_Ntot)
    v <- (logit_v)
    v_a <- exp(f_va)
    v_b <- exp(f_vb)
    alpha <- exp(f_alpha)
    # a <- exp(f_a)
    # b <- exp(f_b)

    ## Negative log-likelihood
    nll <- 0

    nll <- nll - sum(dnorm(v, 0, 1, log = TRUE));       # Assign N(0,1) distribution u
    v_u <-  pnorm(v, 0, 1);  # Uniformly distributed variables (on [0,1])
    v_w <- qbeta(v_u, shape1 = v_a, shape2 = v_b);

    nll <- nll - sum(dnorm(log_delta, 0, 1,log = TRUE));       # Assign N(0,1) distribution u
    delta_u <-  pnorm(log_delta,0, 1);  # Uniformly distributed variables (on [0,1])
    delta_w <- qgamma(delta_u, shape = alpha, scale = 1.0);


    pent <- rep(0,s)
    lambda <- rep(0,s)
    temp <- rep(0,s)
    tau <- rep(0,s)
    psi <- rep(0,s)
    multP <- rep(0,s)

    # Forward algorithm replacement for 2-state model
    pent <- delta_w / sum(delta_w)  # Entry probabilities

    alive <- rep(0, s)
    dead <- rep(0, s)
    temp <- rep(0, s)

    Nsuper <- exp(log_Ntot)  # or (Ntot) depending on your parameterization

    alive[1] <- Nsuper * pent[1]
    dead[1]  <- 0

    for (t in 2:s) {
      survived   <- alive[t - 1] * phi * (1 - p)
      recaptured <- alive[t - 1] * phi * p

      alive[t] <- survived + Nsuper * pent[t]
      dead[t]  <- dead[t - 1] + recaptured
    }

    psi_i <- numeric(s)


    multP_raw <- numeric(s)
    for (t in 1:s) {
      total <- 0
      for (e in 1:(t-1)) {
        alive <- 1
        for (k in (e+1):(t-1)) {
          alive <- alive * phi * (1 - p)
        }
        pr_detected_now <- alive * phi * p
        total <- total + pent[e] * pr_detected_now
      }
      multP_raw[t] <- total
    }

    psiPtot <- sum(multP_raw)  # Make sure this matches your earlier psiPtot

    multP <- multP_raw / psiPtot

    uTot <- sum(u)

    for (i in 1:length(t_k)) {
      entry <- t_k[i]
      recapture <- r_k[i]
      is_tagged <- tag[i]
      count <- n[i]

      p_marked <- v_w[entry]

      if (is_tagged) {
        nll <- nll - dbinom(count, size = count, prob = p_marked, log = TRUE)

        # Marginal recapture probability (forward from entry)
        recapture_prob <- 0
        tmp_alive <- 1
        for (t in (entry + 1):s) {
          tmp_alive <- tmp_alive * phi * (1 - p)
          recapture_prob <- recapture_prob + tmp_alive * p
        }

        if (!is.na(recapture)) {
          nll <- nll - dbinom(count, size = count, prob = recapture_prob, log = TRUE)
        } else {
          nll <- nll - dbinom(0, size = count, prob = recapture_prob, log = TRUE)
        }

      } else {
        nll <- nll - dbinom(count, size = count, prob = 1 - p_marked, log = TRUE)
      }
    }

    nll <- nll - dpois(uTot, Nsuper * psiPtot, log = TRUE)
    nll <- nll - dmultinom(u, prob = multP, log = TRUE)

    out <- list(
      Nsuper = Nsuper,
      v_w = v_w,
      phi = phi,
      lambda = lambda,
      tau = tau,
      p = p,
      psi = psi,
      log_delta = log_delta,
      pent = pent,
      psiPtot = psiPtot,
      multP = multP,
      delta_w = delta_w,
      temp = temp
      # pent_tmp = pent_tmp,
      # psi_tmp = psi_tmp
      # delta_w = exp(log_delta)
    )
    log_p <- log(p)
    log_phi <- log(phi)
    log_psiPtot <- log(psiPtot)
  RTMB::REPORT(out)
  ADREPORT(Nsuper)
  ADREPORT(log_p)
  ADREPORT(log_phi)
  ADREPORT(log_psiPtot)
  return(nll)

}

obj <- RTMB::MakeADFun(model,
                       parameters,
                       map = list(f_va = as.factor(NA)
                                  , f_vb = as.factor(NA)
                                  , f_alpha = as.factor(NA)
                                  ),
                       random = c( "logit_v", "log_delta"),
                       silent = FALSE)

# Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
tmb_rep <- obj$report()
tmb_sdr <- sdreport(obj)
sdr_est <- as.list(tmb_sdr, "Estimate", report=TRUE)
sdr_se <- as.list(tmb_sdr, "Std. Error", report=TRUE)
