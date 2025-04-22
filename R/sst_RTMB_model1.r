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
  summarise(n = n()) #%>%
  # filter((t_l<=r_l) | is.na(r_l))



if(is.null(dim(MR_data$R)[2])==TRUE){
  MR_data[lapply(MR_data,class)=="matrix"]<-lapply(MR_data[lapply(MR_data,class)=="matrix"],as.vector)
  MR_data$alpha<-rep(1,MR_data$s)
}

u <- rep(0,max(d$t_k))
u_i <- aggregate(list(u = d$n),by = list(s = d$t_k),sum)
u[u_i$s] <- u_i$u

sp_count <- d %>%
  group_by(t_l) %>%
  summarise(n = sum(n))

data <- list(
  s = max(d$t_k),                    # Number of periods
  t_k = d$t_k,               # First detection
  r_k = d$r_k,               # Recapture time (NA if none)
  t_l = d$t_l,               # First detection
  r_l = d$r_l,               # Recapture time (NA if none)
  tag = as.integer(d$tag),   # 1 = tagged, 0 = untagged
  n = d$n,                   # Number of fish in each group
  u = u,          # First detections per period (for multinomial)
  uTot = sum(d$n),           # Total number of detections
  mod = 1, #3 forward matrix (1D), 4 forward matrix (multi D)
  sp_pr = as.vector(sp_count$n/sum(sp_count$n))
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

    ## Negative log-likelihood
    nll <- 0

    #Temporal marking probability
    nll <- nll - sum(dnorm(v, 0, 1, log = TRUE));       # Assign N(0,1) distribution u
    v_u <-  pnorm(v, 0, 1);  # Uniformly distributed variables (on [0,1])
    v_w <- qbeta(v_u, shape1 = v_a, shape2 = v_b);


    #Entry probability
    nll <- nll - sum(dnorm(log_delta, 0, 1,log = TRUE));       # Assign N(0,1) distribution u
    delta_u <-  pnorm(log_delta,0, 1);  # Uniformly distributed variables (on [0,1])
    delta_w <- qgamma(delta_u, shape = alpha, scale = 1.0);


    #All the intermediate variables
    pent <- rep(0,s)
    lambda <- rep(0,s)
    temp <- rep(0,s)
    tau <- rep(0,s)
    psi <- rep(0,s)
    multP <- rep(0,s)

    for (i in 1:length(t_k)) {
      entry <- t_k[i]
      recapture <- r_k[i]
      is_tagged <- tag[i]
      count <- n[i]

      p_marked <- v_w[entry]

      lambda[s] <- 0
      for (j in (s-1):1){
        lambda[j] <- phi*(p+(1-p)*lambda[j+1])
      }

      prob <- 0
      alive <- 1
      if (entry < s) {
        for (t in (entry + 1):s) {
          prob <- prob + alive * phi * p
          alive <- alive * phi * (1 - p)
        }
      } else {
        prob <- 0
      }

      if (is_tagged) {
        # These fish were all tagged at first detection, so:
        nll <- nll - dbinom(count, size = count, prob = p_marked, log = TRUE)

        if (!is.na(recapture)) {
          # lambda_fwd <- alpha[2]
          nll <- nll - dbinom(count, size = count, prob = prob, log = TRUE)
        } else {
          # lambda_fwd <- alpha[2]
          nll <- nll - dbinom(0, size = count, prob = prob, log = TRUE)
        }
      } else {
        # These fish were detected but not tagged — they all failed to be tagged
        nll <- nll - dbinom(count - 0, size = count, prob = 1 - p_marked, log = TRUE)
      }
    }

    # if(mod %in% c(1,2)){
    pent <- delta_w/sum(delta_w) #scaled probability of entry
    psi[1] <- pent[1]
    for (i in 1:s){
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

    nll <- nll -
      dpois(uTot, Nsuper * psiPtot , log = TRUE) -
      dmultinom(u, prob = multP ,    log = TRUE)

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
rep <- obj$report()
print(rep$out$Nsuper)

lambda_df <- tibble(
  time = 1:length(rep$out$lambda),
  lambda = rep$out$lambda,
  detections = data$u
)
ggplot(lambda_df, aes(x = time, y = lambda)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(aes(size = detections), color = "darkred", alpha = 0.7) +
  scale_size_continuous(name = "Detections") +
  labs(title = "Recapture Probability (λ) by Time",
       x = "Time (Week or Period)",
       y = "λ (Pr[Recapture | Tagging])") +
  theme_minimal()

df_plot <- tibble(
  time = 1:length(data$u),
  observed = data$u,
  expected = rep$out$multP * sum(data$u)
)

ggplot(df_plot, aes(x = time)) +
  geom_bar(aes(y = observed), stat = "identity", fill = "grey70", alpha = 0.6) +
  geom_line(aes(y = expected), color = "steelblue", size = 1.2) +
  labs(title = "Observed vs. Predicted Weekly Detections",
       y = "Counts", x = "Time (Week or Period)") +
  theme_minimal()

psi_df <- tibble(
  time = 1:length(rep$out$psi),
  psi = rep$out$psi,
  temp = rep$out$temp
)

ggplot(psi_df, aes(x = time)) +
  geom_line(aes(y = psi, color = "psi")) +
  geom_line(aes(y = temp, color = "temp")) +
  labs(title = "Survival and Detection (psi, temp)",
       y = "Value", color = "") +
  theme_minimal()

tibble(
  time = 1:length(rep$out$tau),
  tau = rep$out$tau
) %>%
  ggplot(aes(x = time, y = tau)) +
  geom_line(color = "darkblue", size = 1.2) +
  labs(title = "Conditional Recapture Probability (tau)",
       y = "Probability", x = "Time") +
  theme_minimal()
