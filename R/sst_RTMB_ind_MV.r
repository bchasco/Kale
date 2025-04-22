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
  summarise(n = n()) #%>%
  # filter((t_l<=r_l) | is.na(r_l))



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
  logit_v = rep(0, data$s),   # logit marking rate v[i]
  theta_par = rep(0,15-5)
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

      alive <- 1
      prob <- 0

      lambda[s] <- 0
      for (j in (s-1):1){
        lambda[j] <- phi*(p+(1-p)*lambda[j+1])
      }

      for (t in (entry+1):s) {
        prob <- prob + alive * phi * p               # detection happens *this* step
        alive <- alive * phi * (1 - p)               # update for next step
      }

      if(mod==1){
        prob <- prob
        if(entry==s){
          prob <- 0
        }
      }
      if(mod==2){
        prob <- lambda[entry]
      }
      if(mod==3){
        alpha <- matrix(c(1, 0, 0), nrow = 1)  # Starts in alive state

        S <- matrix(c(
          phi,   0,     1 - phi,
          0,     1,     0,
          0,     0,     1
        ), nrow = 3, byrow = TRUE)

        D <- matrix(c(
          1 - p,  p,   0,
          0,     1,   0,
          0,     0,   1
        ), nrow = 3, byrow = TRUE)

        T <- S %*% D

        for (t in (entry + 1):s) {
          alpha <- alpha %*% T
        }

        prob <- if (entry == s) 0 else as.numeric(alpha[1, 2])
      }
      if(mod==4){
        n_states <- 7
        alpha <- matrix(0, nrow = 1, ncol = n_states)
        alpha[ t_l[i] ] <- 1  # starts alive at its tagging location

        S <- diag(7)
        for (ii in 1:max(na.omit(r_l))) {
          for (j in 1:max(na.omit(r_l))) {
            if(ii == j){
              S[ii, j] <- phi  # survival
            }
          }
          S[ii, 7] <- 1 - phi  # dies without being detected
        }


        D <- diag(7)
        for (l in 1:5) {
          D[l, 6] <- p   # probability of detection from location l
          D[l, l] <- 1 - p
        }

        T <- S %*% D
        alpha[t_l[i]] <- 1  # starts at its tagging location
        for (t in (t_k[i] + 1):s) {
          alpha <- alpha %*% T
        }
        prob <- if (entry == s) 0 else as.numeric(alpha[1,6])
      }
      # if(mod==5){
      #   n_states <- 7
      #   alpha <- matrix(0, nrow = 1, ncol = n_states)
      #   alpha[ t_l[i] ] <- 1  # starts alive at its tagging location
      #
      #   S <- diag(7)
      #   for (ii in 1:max(na.omit(r_l))) {
      #     for (j in 1:max(na.omit(r_l))) {
      #       if(ii == j){
      #         S[ii, j] <- phi  # survival
      #       }
      #     }
      #     S[ii, 7] <- 1 - phi  # dies without being detected
      #   }
      #
      #
      #   D <- diag(7)
      #   for (l in 1:5) {
      #     D[l, 6] <- p   # probability of detection from location l
      #     D[l, l] <- 1 - p
      #   }
      #
      #   T <- S %*% D
      #   alpha[t_l[i]] <- 1  # starts at its tagging location
      #   for (t in (t_k[i] + 1):s) {
      #     alpha <- alpha %*% T
      #   }
      #   prob <- if (entry == s) 0 else as.numeric(alpha[1,6])
      # }
      #
      # if(mod==5){
      #
      #   ## raw movement “scores”
      #   ## diagonal + strictly‑upper cells: 5 + 4 + 3 + 2 + 1 = 15
      #   # theta_par <- rep(0, 15)      # put this in your parameters list
      #
      #   # vec2upper <- function(v) {
      #   #   M <- matrix(0, 5, 5)
      #   #   M[upper.tri(M, diag = TRUE)] <- v
      #   #   M
      #   # }
      #
      #   theta_raw <- diag(5)           # 5×5 with zeros below diag
      #
      #   ## stabilise row‑wise to avoid overflow / underflow
      #   icnt <- 1
      #   for(ii in 1:5){
      #     for(jj in 1:5){
      #       if(jj>ii){
      #         theta_raw[ii,jj] <- exp(theta_par[icnt])
      #         icnt <- icnt + 1
      #       }
      #     }
      #   }
      #   for (r in 1:5) {
      #     theta_raw[r, ] <- theta_raw[r, ]/sum(theta_raw[r, ])
      #   }
      #
      #   # Psi  <- exp(theta_raw)
      #   # Psi[lower.tri(Psi)] <- 0                    # keep them exactly zero
      #   Psi  <- theta_raw
      #
      #   n_states <- 7
      #   S <- matrix(0, 7, 7)                        # 5 live + detected + dead
      #   for (l in 1:5) {
      #     S[l, l:5] <-  phi * Psi[l, l:5]            # survive & (possibly) move up
      #     S[l,     7] <- 1 - phi                     # die unseen
      #   }
      #
      #   D <- diag(7)
      #   for (l in 1:5) {
      #     D[l, 6] <- p                             # detected at ℓ
      #     D[l, l] <- 1 - p
      #   }
      #
      #   T <- S %*% D                                # one‑step kernel
      #
      #   alpha[t_l[i]] <- 1  # starts at its tagging location
      #   for (t in (t_k[i] + 1):s) {
      #     alpha <- alpha %*% T
      #   }
      #
      #
      #   if (entry == s) prob <- 0
      #   if(!is.na(r_l[i])){
      #     prob <- as.numeric(alpha[1,  r_l[i]])
      #   }
      #   if(is.na(r_l[i])){
      #     prob <- as.numeric(alpha[1,  6])
      #   }
      #   # if()
      #   # prob <- if (entry == s) 0 else as.numeric(alpha[1,6])
      # }

      # prob <- if (entry == s) 0 else as.numeric(alpha[1, 2])
      if(mod==1){
        prob <- prob
        if(entry==s){
          prob <- 0
        }
      }
      if(mod==2){
        prob <- lambda[entry]
      }
      if(mod==3){
        prob <- prob
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
      print("Psi mod 1:3")
      print(psi)
      print("temp mod 1:3")
      print(temp)
      print("psiPtot 1")
      print(sum(temp))
      # }

    if(mod == 3){
      pent <- delta_w / sum(delta_w)  # entry probabilities (per time period)
      tmp <- matrix(c(pent[1], 0, 0), nrow = 1)
      psi[1] <- tmp[1,1]  # alive and available (not yet detected)
      for (i in 1:s) {
        # Add new entries before advancing
        tmp <- tmp + matrix(c(pent[i] * (phi - 1)/log(phi), 0, 0), nrow = 1)

        psi[i] <- tmp[1,1]  # alive and available (not yet detected)

        # Propagate to next step
        temp[i] <- psi[i] * p
        tmp <- tmp %*% T
      }
    }

    #Proportion of the total population that is detected. Sum of the joint psi, p, and pent
    if(mod==4){

      print("**************************")
      n_states <- 7                # 5 live locations + detected + dead

      ## build the transition once --------------------------
      S <- diag(n_states)
      for (l in 1:5) S[l, 7] <- 1 - phi          # die unseen
      diag(S)[1:5] <- phi                       # survive in place

      D <- diag(n_states)
      for (l in 1:5) {                          # l = live location
        D[l, 6] <- p                            # detected at l
        D[l, l] <- 1 - p                        # not detected
      }
      T <- S %*% D                              # one‑step kernel

      ## temporary stores -----------------------------------
      psi2  <- matrix(0, s, 5)
      temp2 <- matrix(0, s, 5)

      ## loop over locations j (= 1…5) ----------------------
      ## loop over tagging locations j = 1…5
      for (j in 1:5) {

        alpha <- matrix(0, 1, n_states)
        alpha[1, j] <- pent[1] * sp_pr[j]          # entries that arrive in week 1

        for (t in 1:s) {

          ## add *new* entries that show up in week t
          alpha[1, j] <- alpha[1, j] +
            pent[t] * sp_pr[j] * (phi - 1) / log(phi)

          ## record “alive & still undetected” and the detection term
          psi2 [t, j] <- sum(alpha[1, 1:5])
          temp2[t, j] <- psi2[t, j] * p

          ## advance one survey interval
          alpha <- alpha %*% T
        }
      }
      pent2 <- pent%o%sp_pr
      RTMB::REPORT(pent2)
      print("Second psi")
      print(rowSums(psi2))
      RTMB::REPORT(psi2)
      print("Second temp")
      print(rowSums(temp2))
      temp2 <- rowSums(temp2)
      psi2 <- rowSums(psi2)
      multP2 <- temp2/sum(temp2[1:s])
      psiPtot2 <- sum(temp2)
      print("psiPtot second")
      print(sum(temp2))
      print("multP second")
      print(multP2)


    }
      if(mod==5){

        ## raw movement “scores”
        ## diagonal + strictly‑upper cells: 5 + 4 + 3 + 2 + 1 = 15
        # theta_par <- rep(0, 15)      # put this in your parameters list

        theta_raw <- diag(5)           # 5×5 with zeros below diag

        ## stabilise row‑wise to avoid overflow / underflow
        icnt <- 1
        for(ii in 1:5){
          for(jj in 1:5){
            if(jj>ii){
              theta_raw[ii,jj] <- exp(theta_par[icnt])
              icnt <- icnt + 1
            }
          }
        }
        for (r in 1:5) {
          theta_raw[r, ] <- theta_raw[r, ]/sum(theta_raw[r, ])
        }

        # Psi  <- exp(theta_raw)
        # Psi[lower.tri(Psi)] <- 0                    # keep them exactly zero
        Psi  <- theta_raw

        n_states <- 7
        S <- diag(7)                        # 5 live + detected + dead
        for (l in 1:5) {
          S[l, l:5] <-  phi * Psi[l, l:5]            # survive & (possibly) move up
          S[l,     7] <- 1 - phi                     # die unseen
        }

        D <- diag(7)
        for (l in 1:5) {
          D[l, 6] <- p                             # detected at ℓ
          D[l, l] <- 1 - p
        }

        T <- S %*% D                                # one‑step kernel

        alpha[t_l[i]] <- 1  # starts at its tagging location
        for (t in (t_k[i] + 1):s) {
          alpha <- alpha %*% T
        }

        RTMB::REPORT(theta_raw)
        pent <- delta_w / sum(delta_w)  # entry probabilities (per time period)
        # pent2 <- pent %o% sp_pr
        # psi <- numeric(s)
        tmp <- matrix(c(pent[1], 0,0,0,0, 0, 0), nrow = 1)
        psi2 <- matrix(0,length(delta_w),5)
        temp2 <- matrix(0,length(delta_w),5)

        n_states <- 7
        alpha <- matrix(0, nrow = 1, ncol = n_states)
        alpha[ t_l[i] ] <- 1  # starts alive at its tagging location


        for(j in 1:5){
          tmp <- matrix(c(pent[1]*sp_pr[j], 0,0,0,0,0, 0), nrow = 1)
          psi2[1,j] <- tmp[1,1]  # alive and available (not yet detected)
          for (i in 1:s) {
            # Add new entries before advancing
            tmp <- tmp + matrix(c(pent[i]*sp_pr[j] * (phi - 1)/log(phi), 0, 0,0,0,0,0), nrow = 1)

            psi2[i,j] <- tmp[1,1]  # alive and available (not yet detected)

            # Propagate to next step
            temp2[i,j] <- psi2[i,j] * p
            tmp <- tmp %*% T
          }
        }
        pent2 <- pent%o%sp_pr
        RTMB::REPORT(pent2)
        print("Second psi")
        print(rowSums(psi2))
        RTMB::REPORT(psi2)
        print("Second temp")
        print(rowSums(temp2))
        temp2 <- rowSums(temp2)
        psi2 <- rowSums(psi2)
        multP2 <- temp2/sum(temp2[1:s])
        psiPtot2 <- sum(temp2)
        print("psiPtot second")
        print(sum(temp2))
        print("multP second")
        print(multP2)
        RTMB::REPORT(Psi)

        nll <- nll - sum(dnorm(theta_par, 0, 2, log = TRUE))
      }
      if(mod<4){
      #probability of entry
      print('model < 4')
      multP <- temp/sum(temp[1:s])
      psiPtot <- sum(temp)
      print("multP first")
      print(multP)
    }
    if(mod == 4){
      print("model 4")
      multP <- temp2/sum(temp2[1:s])
      psiPtot <- sum(temp2)
      print("multP first")
      print(multP)
    }
      if(mod == 5){
        print("model 4")
        multP <- temp2/sum(temp2[1:s])
        psiPtot <- sum(temp2)
        print("multP first")
        print(multP)
      }

    uTot <- sum(u)
    Nsuper <- (Ntot)

    if (mod == 4 | mod == 5) {
      nll <- nll -
        dpois(uTot, Nsuper * psiPtot2, log = TRUE) -
        dmultinom(u, prob = multP2,    log = TRUE)
    } else {
      nll <- nll -
        dpois(uTot, Nsuper * psiPtot , log = TRUE) -
        dmultinom(u, prob = multP ,    log = TRUE)
    }

    lambda_scalar <- numeric(s)
    lambda_matrix <- numeric(s)
    lambda_matrix2 <- numeric(s)

    for (entry in 1:(s-1)) {
      # Scalar method
      alive <- 1
      prob <- 0
      for (t in (entry + 1):s) {
        prob <- prob + alive * phi * p
        alive <- alive * phi * (1 - p)
      }
      lambda_scalar[entry] <- prob

      # Matrix method
      alpha <- matrix(c(1, 0), nrow = 1)
      T <- matrix(c(phi * (1 - p), phi * p,
                    0           , 1), nrow = 2, byrow = TRUE)
      for (t in (entry + 1):s) {
        alpha <- alpha %*% T
      }
      lambda_matrix[entry] <- alpha[1, 2]

  }

    lambda_matrix3 <- numeric(s)
    for (entry in 1:(s-1)) {
      alpha <- matrix(c(1, 0, 0), nrow = 1)
      S <- matrix(c(
        phi,   0,     1 - phi,
        0,     1,     0,
        0,     0,     1
      ), nrow = 3, byrow = TRUE)

      D <- matrix(c(
        1 - p,  p,   0,
        0,     1,   0,
        0,     0,   1
      ), nrow = 3, byrow = TRUE)

      T <- S %*% D
      for (t in (entry + 1):s) {
        alpha <- alpha %*% T
      }

      lambda_matrix3[entry] <- alpha[1, 2]
    }

    lambda_matrix4 <- numeric(s)
    if(mod==4){
      n_states <- 7
      for (entry in 1:(s-1)) {
        alpha <- matrix(0, nrow = 1, ncol = n_states)
        alpha[ 2 ] <- 1  # starts alive at its tagging location

        S <- diag(7)
        for (i in 1:max(na.omit(r_l))) {
          for (j in 1:max(na.omit(r_l))) {
            if(i == j){
              S[i, j] <- phi  # survival
            }
          }
          S[i, 7] <- 1 - phi  # dies without being detected
        }


        D <- diag(7)
        for (l in 1:5) {
          D[l, 6] <- p   # probability of detection from location l
          D[l, l] <- 1 - p
        }

        T <- S %*% D
        for (t in (entry + 1):s) {
          alpha <- alpha %*% T
        }
        lambda_matrix4[entry] <- alpha[1, 6]
      }
    }

    if(mod==5){

        ## raw movement “scores”
        ## diagonal + strictly‑upper cells: 5 + 4 + 3 + 2 + 1 = 15
        # theta_par <- rep(0, 15)      # put this in your parameters list

        theta_raw <- diag(5)           # 5×5 with zeros below diag

        ## stabilise row‑wise to avoid overflow / underflow
        icnt <- 1
        for(ii in 1:5){
          for(jj in 1:5){
            if(jj>ii){
              theta_raw[ii,jj] <- exp(theta_par[icnt])
              icnt <- icnt + 1
            }
          }
        }
        for (r in 1:5) {
          theta_raw[r, ] <- theta_raw[r, ]/sum(theta_raw[r, ])
        }

        # Psi  <- exp(theta_raw)
        # Psi[lower.tri(Psi)] <- 0                    # keep them exactly zero
        Psi  <- theta_raw

        n_states <- 7
        S <- diag(7)                        # 5 live + detected + dead
        for (l in 1:5) {
          S[l, l:5] <-  phi * Psi[l, l:5]            # survive & (possibly) move up
          S[l,     7] <- 1 - phi                     # die unseen
        }

        D <- diag(7)
        for (l in 1:5) {
          D[l, 6] <- p                             # detected at ℓ
          D[l, l] <- 1 - p
        }
    }
    print("scalar")
    print(lambda_scalar)
    print("matrix")
    print(lambda_matrix)
    print("matrix3")
    print(lambda_matrix3)
    print("matrix4")
    print(lambda_matrix4)

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
      temp = temp,
      S = S,
      T = T,
      D = D
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
# tmb_sdr <- sdreport(obj)
# sdr_est <- as.list(tmb_sdr, "Estimate", report=TRUE)
# sdr_se <- as.list(tmb_sdr, "Std. Error", report=TRUE)

print(rep$out$Nsuper)
print(rep$out$phi)
print(rep$out$p)
plot(data$u/sum(data$u))
lines(rep$out$pent)
# lines(1:length(rep$out$v_w)+0.1, rowSums(rep$pent2))
