library(fs)
library(tibble)
library(knitr)
library(dplyr)
library(tidyr)
library(readr)

# Define the relative directory (from project root)
directory_path <- "Data/nfl"

# Get the list of files (returning only filenames)
file_list <- dir_ls(directory_path, type = "file") |>
  basename()

# Function to read each file
read_file_contents <- function(file_path) {
  tryCatch(
    {
      df <- read_csv(paste0("data/nfl/",file_path), col_types = cols(default = "c"), show_col_types = FALSE) |>
        mutate(source_file = basename(file_path)) # Add column for source file

    },
    error = function(e) {
      message("Skipping file: ", file_path, " - Error: ", e$message)
      return(NULL)
    }
  )
}

# Read and combine all files
combined_data <- file_list  %>%
  lapply(read_file_contents)
cols <- c("Return_Yr","SPECIES","Run","TagDate","TagReach",
          "Tag1","Tag2","RecapDate","RecapReach")

result <- do.call(rbind, lapply(combined_data, function(df) df[, cols, drop = FALSE])) %>%
  mutate(t_wk = lubridate::week(lubridate::dmy(TagDate)),
         r_wk = ifelse(is.na(RecapDate),NA,lubridate::week(lubridate::dmy(RecapDate))),
         y_i = Return_Yr) %>%
  mutate(r_wk = ifelse(as.numeric(as.character(r_wk))<20,r_wk+ 52,r_wk),
         t_wk = ifelse(as.numeric(as.character(t_wk))<20,t_wk+ 52,t_wk)) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1)

# Create lookup tables separately for TagReach and RecapReach
lk_up <- data.frame(
  TagReach = sort(unique(result$TagReach)),
  t_l = seq_along(unique(result$TagReach))
)

lk_up_r <- data.frame(
  TagReach = sort(unique(result$TagReach)),
  r_l = seq_along(unique(result$TagReach))
)

# Perform left joins to assign t_l and r_l separately
d <- result %>%
  filter(t_k<=22,
         Return_Yr %in% c(2014)) %>%
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




# d <- read.csv("data/simpleData2.csv") %>%
#   mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
#          r_wk = lubridate::week(lubridate::mdy(RecapDate)),
#          t_yr = lubridate::year(lubridate::mdy(TagDate)),
#          r_yr = lubridate::year(lubridate::mdy(RecapDate))) %>%
#   filter(t_yr == 2023) %>%
#   filter(t_wk > 10) %>%
#   mutate(t_k = t_wk - min(t_wk) + 1,
#          r_k = r_wk - min(t_wk) + 1) %>%
#   mutate(t_l = TagState,
#          r_l = RecapState) %>%
#   mutate(y_i = as.integer(factor(Return_Yr))) %>%
#   filter(is.na(r_wk) | r_k>0) %>%
#   mutate(tag = ifelse(Tag1=="",FALSE,TRUE)) %>%
#   group_by(y_i,t_k,r_k,t_l, r_l, tag) %>%
#   summarise(n = n())

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
  surv_par = matrix(-2,data$n_y,max_l * (max_l + 1) / 2 - max_l),      # Logit persistence probability
  mu_detect = 0,
  detection_par = matrix(0,data$n_y,max_l-1),        # Logit detection probability
  # detection_re = array(0,c(data$n_y,5,max(data$t_k))),      # Logit tagging probability
  taggingRate_par = matrix(0,data$n_y,max_l-1),      # Logit tagging probability
  taggingRate_re = array(0,c(data$n_y,max_l-1,max(data$t_k))),      # Logit tagging probability
  par_PopTotal = rep(0,data$n_y),   # Log expected carcass births
  B_time_sig = 0,
  B_loc_sig = 0,
  taggingRate_sig = 0,
  detection_sig = 0,
  B_time = matrix(0, data$n_y, max(data$t_k)),  #Carcasses by week
  B_loc = array(0,c(data$n_y, max_l-1, max(data$t_k)))  #Distribution of carcasses by location
)

# Define RTMB Model
rtmb_model <- function(parms){

    RTMB::getAll(data, parms)

    nll <- 0

    max_l <- max(na.omit(r_l)) + 1

    print(max_l)

    nll_detect <- array(1,c(n_y,max_l,max(t_k)))
    nll_tag <- array(0,c(n_y,max_l,max(t_k)))
    nll_CJS <- rep(0,length(t_k))

    #Survival
    surv <- array(0,c(n_y,max_l,max_l))
    f_surv <- array(0,c(n_y,max_l,max_l))
    #Survival
    for(y in unique(y_i)){
      ii <- 1
      for(i in 1:(max_l-1)){
        for(j in (i+1):max_l){
          surv[y,i,j] <- exp(surv_par[y,ii])
          f_surv[y,i,j] <- exp(surv_par[y,ii])
          ii <- ii + 1
        }
        surv[y,i,i] <- -sum(surv[y,i,])
        f_surv[y,i,i] <- -sum(f_surv[y,i,])
      }
      surv[y,,] <- as.matrix(Matrix::expm(surv[y,,]))

      # surv2[y,,] <- t(surv[y,,])
    }

    # Compute probability matrix
    p <- array(0,c(n_y,max(t_k),max_l,max_l))
    for(y in unique(y_i)){
      for(t in 1:max(t_k)){
        diag(p[y,t,,]) <- c(plogis(mu_detect + detection_par[y,]),0)
        p[y,t,,max_l] <- 1 - diag(p[y,t,,])
        p[y,t,max_l,max_l] <- 1
      }
    }

    # Derived variables
    N <- array(0, c(n_y,max_l,max(t_k)))      # Total available carcasses over time
    TotalCarcasses_t <- array(0,c(n_y,max_l,max(t_k)))      # Observed Total carcasses over time
    TaggedCarcasses_t <- array(0,c(n_y,max_l,max(t_k)))      # Observed Tagged carcasses over time
    E_TotalCarcasses_t <- array(0,c(n_y,max_l,max(t_k)))      # Expected Total carcasses over time
    E_TaggedCarcasses_t <- array(0,c(n_y,max_l,max(t_k)))      # Expected Tagged carcasses over time

    #Conditional probability of loc for each carcass
    B_loc_dist <- array(0,c(n_y,max_l,max(t_k)))
    for(y in unique(y_i)){
      B_loc_dist[y,,1] <- exp(c(B_loc[y,,1],0))/sum(exp(c(B_loc[y,,1],0)))
      # Initialize spatial distribution for the first time step
      N[y,,1] <- exp(B_time[y,1] + par_PopTotal[y]) * B_loc_dist[y,,1]
      # Loop through time steps to compute derived variables
      N[y,,t] = t(Matrix::expm(f_surv[y,,] * 0.5)) %*% (exp(B_time[y,t] + par_PopTotal[y]) * B_loc_dist[y,,t]) +
        t(surv[y,,]) %*% N[y,,t-1]
    }

    B_ts <- B_loc_dist * 0
    B_total <- rep(0,n_y)
    for(y in unique(y_i)){
      B_ts[y,,] <- log(t(B[y,] * t(B_loc_dist[y,,])))
      B_total[y] <- log(sum(B[y,]))
    }

    print(B_total)
    print("test")
    B <- exp(B_time + par_PopTotal)
    nll_B <- rep(0,max(t_k))

    for(y in unique(y_i)){
      for (t in unique(t_k)) {
        # Detection Process (Binomial likelihood)
        for(loc in 1:(max_l-1)){
          TaggedCarcasses_t[y,loc,t] <- sum(na.omit(n[y_i == y & t_k == t & t_l == loc & tag == 1]))  # Tagged carcasses at time t
          TotalCarcasses_t[y,loc,t] <- sum(na.omit(n[y_i == y & t_k == t & t_l == loc]))  # Total detected carcasses at time t
        }

        #Carcass detection probability likelihood
        E_TotalCarcasses_t[y,1:(max_l-1),t] <- N[y,1:(max_l-1),t]*plogis(mu_detect + detection_par[y,])
        nll_detect[y,1:(max_l-1),t] <- RTMB::dpois(TotalCarcasses_t[y,1:(max_l-1),t]+1,
                                         E_TotalCarcasses_t[y,1:(max_l-1),t],
                                         log = TRUE)

        #Carcass tagging rate likelihood
        for(loc in 1:(max_l-1)){
          E_TaggedCarcasses_t[y,loc,t] <- E_TotalCarcasses_t[y,loc,t] * plogis(taggingRate_par[y,loc] + taggingRate_re[y,loc,t])
          if(TotalCarcasses_t[y,loc,t]>0)
            nll_tag[y,loc,t] <- RTMB::dpois(TaggedCarcasses_t[y,loc,t] + 1,
                                            E_TaggedCarcasses_t[y,loc,t],
                                          log = TRUE)
        }
      }
    }



    for(i in 1:length(t_l)){
      print(i)
      # if(tag[i]==TRUE){
        m <- matrix(0,max_l,max_l)
        diag(m) <- 1
        last_k <- ifelse(is.na(r_k[i]),max(na.omit(r_k)),r_k[i])

        if(is.na(r_k[i])){#never recapped
          print("never")
          if((t_k[i]+1)<last_k){
            for(k in (t_k[i]+1):last_k){
              m <- m %*% surv[y_i[i],,] %*% diag(p[y_i[i],k,,max_l])
            }
          }
        }


        if(!is.na(r_k[i])){#Recaptured
          print("recaptured")
          if(r_k[i]>0){
            # print(r_k[i])
            for(k in (t_k[i]+1):r_k[i]){
              if(k<last_k){
                m <- m %*% surv[y_i[i],,] %*% diag(p[y_i[i],k,,max_l]) #non-detection
              }else{
                m <- m %*% surv[y_i[i],,] %*% diag(p[y_i[i],k,,r_l[i]]) #recap
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
      # }
    }

    nll_CJS_total <- sum(log(nll_CJS[tag==TRUE])*n[tag==TRUE])

    nll_B_time <- RTMB::dnorm(B_time, 0, exp(B_time_sig), log = TRUE)
    nll_B_loc <- RTMB::dnorm(B_loc, 0, exp(B_loc_sig), log = TRUE)
    nll_taggingRate_re <- RTMB::dnorm(taggingRate_re, 0, exp(taggingRate_sig), log = TRUE)
    # nll_detection_re <- RTMB::dnorm(detection_re, 0, exp(detection_sig), log = TRUE)
    nll_detection_re <- RTMB::dnorm(detection_par, 0, exp(detection_sig), log = TRUE)

    taggingRate <- taggingRate_re * 0

    for(y in unique(y_i)){
      taggingRate[y,,] <- plogis(taggingRate_re[y,,] + taggingRate_par[y,])
    }

    probability_of_outcome <- nll_CJS

    E <- d %>%
      group_by(t_k, t_l, tag) %>%
      mutate(t_k_total = sum(n))  %>%
      mutate(obs = n)
    E$p <-  probability_of_outcome
    pred2 <- log(E$p * E$t_k_total)


    # REPORT(detection_re)
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
            sum(nll_taggingRate_re)+
              sum(nll_detection_re))
    )
}

obj <- RTMB::MakeADFun(rtmb_model,
                       parameters,
                       random = c("B_time","B_loc","taggingRate_re", "detection_par"),
                       map = list(detection_par = as.factor(rep(NA,length(parameters$detection_par))),
                                  detection_sig = as.factor(NA)
                                  B_loc = as.factor(array(NA,c(data$n_y, max_l-1, max(data$t_k)))),
                                  B_loc_sig = as.factor(NA)),

                       silent = FALSE)

# Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- obj$report()
sdr <- sdreport(obj)
sd.est <- as.list(sdr,"Estimate", report = TRUE)
sd.sd <- as.list(sdr,"Std. Error", report = TRUE)

print(sum(rep$B))


E <- reshape2::melt(sd.est$E_TotalCarcasses_t)
E$sd <- reshape2::melt(sd.sd$E_TotalCarcasses_t)$value
E$obs <- reshape2::melt(rep$TotalCarcasses_t)$value
E %>% filter(Var1 != 6) %>%
  ggplot(aes(x = Var2, y = obs)) +
  geom_point() +
  geom_point(aes(x = Var2, y = value), color = "red", alpha = 0.5) +
  geom_errorbar(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd), color = "red", alpha = 0.5, width = 0.2) +
  facet_grid(Var1~Var3, scales = "free_y") +
  xlab("Location") +
  ylab("Number of total carcasses") +
  theme_classic()

E <- reshape2::melt(sd.est$E_TaggedCarcasses_t)
E$sd <- reshape2::melt(sd.sd$E_TaggedCarcasses_t)$value
E$obs <- reshape2::melt(rep$TaggedCarcasses_t)$value
E %>% filter(Var1 != 6) %>%
  ggplot(aes(x = Var2, y = obs)) +
  geom_point() +
  geom_point(aes(x = Var2, y = value), color = "red", alpha = 0.5) +
  geom_errorbar(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd), color = "red", alpha = 0.5, width = 0.2) +
  facet_grid(Var1~Var3, scales = "free_y") +
  xlab("Location") +
  ylab("Number of tagged carcasses") +
  theme_classic()

E <- d %>%
  group_by(y_i,t_k, t_l, tag) %>%

  mutate(t_k_total = sum(n))  %>%
  mutate(obs = n)
# E$p <-  rep$probability_of_outcome
E$pred <- sd.est$pred2
E$sd <- sd.sd$pred2

E <- E %>%
  filter(tag ==TRUE) %>%
  group_by(y_i,t_k,t_l,r_l) %>%
  summarise(obs = sum(obs),
            pred = sum(exp(pred)),
            sd = sum(sd))

E %>%
  # filter(t_k<=5) %>%
  mutate(r_l = ifelse(is.na(r_l),2,r_l)) %>%
  ggplot(aes(x = r_l, y = obs, color = as.factor(y_i))) +
  geom_point(alpha = 0.8, shape = 2) +
  geom_point(aes(x = r_l, y = (pred), color = as.factor(y_i)), alpha = 0.5) +
  # geom_errorbar(aes(ymin = exp(log(pred) - 1.96 * sd), ymax = exp(log(pred) + 1.96 * sd), color = as.factor(y_i)), alpha = 0.5, width = 0.5) +
  facet_wrap(~t_k, scales = "free_y") +
  xlab("Recapture location") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

E %>%
  filter(t_k>5 & t_k<=10) %>%
  mutate(r_l = ifelse(is.na(r_l),6,r_l)) %>%
  ggplot(aes(x = r_l, y = obs)) +
  geom_point(color = "black", alpha = 0.8) +
  geom_point(aes(x = r_l, y = (pred)), col = "red", alpha = 0.5) +
  # geom_errorbar(aes(ymin = exp(log(pred) - 1.96 * sd), ymax = exp(log(pred) + 1.96 * sd)), col = "red", alpha = 0.5, width = 0.5) +
  facet_grid(t_l~t_k, scales = "free_y") +
  xlab("Recapture location") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

E <- reshape2::melt(sd.est$B_total)
E$sd <- reshape2::melt(sd.sd$B_total)$value
x <- rnorm(1000,sd.est$B_total,sd.sd$B_total)
# hist(exp(x),
#      main = "",
#      ylab = "",
#      xlab = "Total births")
E %>% #filter(Var1 != 6) %>%
  ggplot(aes(x = Var1, y = exp(value))) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(value - 1.96 * sd), ymax = exp(value + 1.96 * sd))) +
  xlab("Year") +
  ylab("Total births") +
  theme_classic()

E <- reshape2::melt(sd.est$B_ts)
E$sd <- reshape2::melt(sd.sd$B_ts)$value
E %>%
  filter(Var2 != 2) %>%
  ggplot(aes(x = Var2, y = exp(value))) +
  geom_point() +
  # geom_errorbar(aes(ymin = exp(value - 1.96 * sd), ymax = exp(value + 1.96 * sd))) +
  facet_grid(Var1~Var3, scales = "free_y") +
  xlab("Location") +
  ylab("Births") +
  theme_classic()

df_long <- reshape2::melt(rep$surv) %>%
  mutate(value = ifelse(value==0,NA,value))
# M <- rep$surv
# # M[6,6] <- NA
# colnames(M) <- c(paste("NLL", 1:5),"'Death'")
# rownames(M) <- c(paste("NLL", 1:5),"'Death'")

ggplot(df_long, aes(x = as.factor(Var3), y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
  facet_wrap(~ Var1, ncol = 2) +
  labs(title = "Trnasition probabilities", fill = "Transition \n probability") +
  theme_minimal() +
  xlab("Initital location") +
  ylab("Recapture location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = df_long %>% filter(!is.na(value)),
            aes(label = round(value, 2)),
            color = "black",
            size = 3)
#
# corrplot::corrplot(
#   M,
#
#   method = "color",       # or "color", "ellipse", etc.
#   type = "upper",          # upper or lower triangle of correlation matrix
#   addCoef.col = "black",   # color of correlation coefficients
#   tl.col = "black",        # color of text labels (variable names)
#   tl.srt = 45,
#   col.lim = c(0,1)
#   # text label rotation
# )

