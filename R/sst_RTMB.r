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

# Create data object for JAGS model
data <-
  list(
    # alpha = 0.5
     LL = sum(collapsed_df$ui)
    , UL = sum(collapsed_df$ui) * 20
    , s = nrow(collapsed_df)
    #, num_strata = dim(n)[2]
    , uTot = sum(collapsed_df$ui)
    , n = collapsed_df$ni
    , R = collapsed_df$Ri
    , u = collapsed_df$ui
    , r = collapsed_df$ri
    ,m = collapsed_df$mi
    , T = collapsed_df$zi + collapsed_df$mi
    , days=c(NA,as.numeric(diff(lubridate::ymd(collapsed_df$dates)))) #number of days between last period and current period...for first period, value is NA
    , good_s_T_mat = good_s_T_mat_2
    , good_s_R_mat = good_s_R_mat_2
    , good_s_T = good_s_T_final
    , good_s_R = good_s_R_final
    , logit_phi_1_mu_prior = logit(exp(log(0.5)/mean(c(NA,as.numeric(diff(lubridate::ymd(collapsed_df$dates)))),na.rm=TRUE)))
  )

if(is.null(dim(MR_data$R)[2])==TRUE){
  MR_data[lapply(MR_data,class)=="matrix"]<-lapply(MR_data[lapply(MR_data,class)=="matrix"],as.vector)
  MR_data$alpha<-rep(1,MR_data$s)
}

# Create inits for uTot
uTot<-MR_data$uTot
inits_func<-function(){
  inits<-list(
    "Ntot" = uTot * runif(1,1,5)
  )
}


## Parameters
parameters <- list(
  f_alpha = log(0.5),
  f_va= qlogis(0.5),
  f_vb= qlogis(0.5),
  logit_phi = 0,        # logit of survival
  logit_p = 0,          # logit of detection
  log_Ntot = 10,        # log total N
  log_delta = rep(-0, data$s),# log of delta (entry intensity)
  logit_v = rep(0, data$s-1)   # logit marking rate v[i]
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

    lambda[s] <- 0
    psi[1] <- delta_w[1]/sum(delta_w)

    # for (i in 1:s){
    #   pent[i] <- delta_w[i]/sum(delta_w) #scaled probability of entry
    # }
    pent <- delta_w/sum(delta_w) #scaled probability of entry
    for (i in (s-1):1){
      lambda[i] <- phi*(p+(1-p)*lambda[i+1])
    }

    psi[1] <- pent[1]
    for (i in 1:s){
      #i is period
      if(i < s){
        psi[i+1] <- psi[i]*(1-p)*phi + pent[i+1] *(phi-1)/log(phi)  #survival and still availabe for capture
      }
      temp[i] <- psi[i]*p #joint survival after entry and then detected
      tau[i] <-p/(p+(1-p)*lambda[i])  #Recapture probabiilty
    }
    multP <- temp/sum(temp[1:s])


    psiPtot <- sum(temp)
    uTot <- sum(u)
    Nsuper <- (Ntot)


    #data
    for (i in 1:(s - 1)) {
      #marks ~ total carcasses
      nll <- nll - dbinom(R[i], size = n[i], prob = v_w[i], log = TRUE)
      #marks ~ marks that were eventually recapped
      nll <- nll - dbinom(r[i], size = R[i], prob = lambda[i], log = TRUE)
    }

    #Data
    for (i in 2:(s - 1)) { #recaps ~ recaps + future recaps
      nll <- nll - dbinom(m[i], size = T[i], prob = tau[i], log = TRUE)
    }


    #priors
    nll <- nll - dbeta(phi, shape1 = v_a, shape2 = v_b, log = TRUE)
    nll <- nll - dbeta(p, shape1 = v_a, shape2 = v_b, log = TRUE)
    # print(nll)

    #observations
    print(psiPtot)
    nll <- nll - dpois(uTot, Nsuper * psiPtot, log = TRUE)
    print(dpois(uTot, Nsuper * psiPtot, log = TRUE))
    nll <- nll - dmultinom(u, prob = multP, log = TRUE)
    print(multP)
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
                       map = list(#f_va = as.factor(NA)
                                  #, f_vb = as.factor(NA)
                                  # , f_alpha = as.factor(NA)
                                  ),
                       random = c( "logit_v", "log_delta"),
                       silent = FALSE)

# Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
tmb_rep <- obj$report()
tmb_sdr <- sdreport(obj)
sdr_est <- as.list(tmb_sdr, "Estimate", report=TRUE)
sdr_se <- as.list(tmb_sdr, "Std. Error", report=TRUE)
