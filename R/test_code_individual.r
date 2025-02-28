library(dplyr)
library(tidyr)

d <- read.csv("data/simpleData2.csv") %>%
  mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
         r_wk = lubridate::week(lubridate::mdy(RecapDate)),
         t_yr = lubridate::year(lubridate::mdy(TagDate)),
         r_yr = lubridate::year(lubridate::mdy(RecapDate))) %>%
  filter(t_yr == 2024) %>%
  filter(t_wk > 10) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  filter(is.na(r_wk) | r_k>0) %>%
  mutate(tag = ifelse(Tag1=="",FALSE,TRUE)) %>%
  group_by(t_k,r_k,tag) %>%
  summarise(n = n())

B_obs <- d %>%
  group_by(t_k) %>%
  summarize(n = sum(n))

T_obs <- d %>%
  filter(tag==TRUE) %>%
  group_by(t_k) %>%
  summarize(n = sum(n))

D_obs <- d %>%
  filter(tag==FALSE) %>%
  group_by(t_k) %>%
  summarize(n = sum(n))

library(RTMB)

# Load Data
data <- list(
  t_k = d$t_k,   # First detection week
  r_k = d$r_k,   # Recapture week (NA if not recaptured)
  tag = as.integer(d$tag), # Whether tagged (1 = tagged, 0 = not)
  n = d$n,        # Number of carcasses in each row
  B_obs = B_obs,
  C_t = T_obs$n,
  D_obs = D_obs
)

parameters = list(
  phi_par = 0,      # Log persistence probability
  p_par = 0,        # Log detection probability
  psi_par = 0,      # Logit tagging probability
  par_lambda = 20,   # Log expected carcass births
  B_sig = 0,
  C_sig = 0,
  R_sig = 0,
  B = rep(2, max(data$t_k))  # Estimated carcass entries (initial guess)
)

# Define RTMB Model
rtmb_model <- function(parms){

    RTMB::getAll(data, parms)

    nll <- 0

    # Define negative log-likelihood components
    nll_birth <- rep(0, length(unique(t_k)))
    nll_detect <- rep(0, length(unique(t_k)))
    nll_tag <- rep(0, length(unique(t_k)))
    nll_recapture <- rep(0, length(unique(t_k)))


    # Derived variables
    N <- rep(0,max(t_k))      # Total carcasses over time
    C_t <- rep(0,max(t_k))      # Total carcasses present at time t
    T_t <- rep(0,max(t_k))      # Fish tagged at time t
    T_available <- rep(0,max(t_k))  # Tagged carcasses over time
    E_R <- matrix(0,max(t_k),max(t_k))  # Expected recaptures
    R_t <- matrix(0,max(t_k),max(t_k))  # Observed recaptures

    # Initialize first time step
    N[1] <- exp(B[1] + par_lambda)
    T_available[1] <- exp(B[1] + par_lambda) * plogis(psi_par);

    # Loop through time steps to compute derived variables
    for (t in 2:max(t_k)) {

      # Total carcasses at time t (new + survivors)
      N[t] = exp(B[t] + par_lambda) + N[t-1] * plogis(phi_par);

      # Tagged carcasses at time t (new tags + surviving tags)
      T_available[t] = sum(n[t_k == t & tag == 1]) + T_available[t-1] * plogis(phi_par);
    }

    #Number of births
    nll_birth <- RTMB::dnorm(B, 0, exp(B_sig), log = TRUE)

    for (t in unique(t_k)) {

      # Detection Process (Binomial likelihood)
      C_t[t] <- sum(n[t_k == t])  # Total detected carcasses at time t
      nll_detect[t] <- RTMB::dpois(C_t[t], N[t]*plogis(p_par), log = TRUE)

      # Tagging Process (Binomial likelihood)
      #Conditional on the observed carcasses, how many are tagged
      T_t[t] <- sum(n[t_k == t & tag == 1])  # Tagged carcasses at time t
      nll_tag[t] <- RTMB::dbinom(T_t[t], C_t[t], plogis(psi_par), log = TRUE)


      # Recapture Process (Poissson )
      if (t > 1 & t< 16) {
        for (r in (t + 1):max(t_k)) {  # Future time steps

          E_R[t, r] <- T_available[t] * plogis(phi_par)^(r - t)
          R_t[t,r] <- sum(na.omit(n[r_k == r & t_k == t & tag == 1]))  # Recaptured tagged carcasses
          # Likelihood for observed recaptures
          if(R_t[t,r]>0){
            nll_recapture[t] <- nll_recapture[t] +  RTMB::dpois(R_t[t,r], E_R[t, r]  * plogis(p_par), log = TRUE)
          }

        }

        # # #All the fish that were never observed
        # tmp_Rt[t] <- sum(na.omit(n[is.na(r_k) & t_k == t & tag == 1]))
        # E_R_tmp[t] <- T_available[t] * prod(1-plogis(phi_par)^(t:16 - (t-1)) * plogis(p_par))
        # nll_recapture[t] <- nll_recapture[t] + RTMB::dpois(tmp_Rt[t], E_R_tmp[t], log = TRUE)
      }
    }

    # Combine all likelihood components
    nll <- -sum(nll_birth, na.rm = TRUE) -
      sum(nll_detect, na.rm = TRUE) -
      sum(nll_tag, na.rm = TRUE) -
      sum(nll_recapture, na.rm = TRUE)

    B_star <- exp(B + par_lambda)
    N_total <- sum(B_star)

    RTMB::REPORT(B)
    RTMB::REPORT(B_star)
    RTMB::REPORT(T_t)
    RTMB::REPORT(tmp_Rt)
    RTMB::REPORT(E_R_tmp)
    RTMB::REPORT(C_t)
    RTMB::REPORT(R_t)
    RTMB::REPORT(E_R)
    RTMB::REPORT(E_Rp)
    RTMB::REPORT(N)
    RTMB::REPORT(T_available)
    RTMB::REPORT(nll_recapture)
    RTMB::REPORT(nll_tag)
    RTMB::REPORT(nll_detect)
    RTMB::REPORT(nll_birth)
    RTMB::ADREPORT(B_star)
    RTMB::ADREPORT(N_total)
    return(nll)
}

obj <- RTMB::MakeADFun(rtmb_model,
                       parameters,
                       random = c("B"),
                       map = list(
                         R_sig = as.factor(NA),
                         C_sig = as.factor(NA)),
                       silent = FALSE)

# Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- obj$report()
sd <- sdreport(obj)

E <- reshape2::melt(rep$E_R)
E$obs <- reshape2::melt(rep$R_t)$value

library(ggplot2)
E %>% ggplot(aes(x = Var1, y = obs)) +
  geom_point() +
  geom_line(aes(x = Var1, y = value)) +
  facet_wrap(~Var2, ncol = 4) +
  theme_classic()


# plot(rep$N)
