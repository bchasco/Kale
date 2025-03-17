year <- c(2021)


  library(fs)
  library(tibble)
  library(knitr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)

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
  cols <- c("Return_Yr","SPECIES","Run","ScaleAge", "Sex","TagDate","TagReach",
            "Tag1","Tag2","RecapDate","RecapReach")

  result <- do.call(rbind, lapply(combined_data, function(df) df[, cols, drop = FALSE])) %>%
    mutate(oa = ifelse(is.na(ScaleAge)|(ScaleAge==9)|(ScaleAge==88),-1,as.numeric(substr(ScaleAge,1,1)) - as.numeric(substr(ScaleAge,2,2)))) %>%
    mutate(Sex = ifelse(is.na(Sex)|(Sex=="Unknown")|(Sex=="Adult"),-1,Sex)) %>%
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
  #
  # Perform left joins to assign t_l and r_l separately
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
    group_by(y_i,t_k,r_k,Sex,TagReach,RecapReach,t_l, r_l, tag) %>%
    summarise(n = n()) %>%
    group_by(y_i, Sex) %>%
    mutate(grp = cur_group_id()-1, #The minus 1 removes the -1 effect
           grp_id = paste(y_i,Sex)) %>%
    ungroup()


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
    grp = d$grp,
    n_y = length(unique(d$y_i))
  )

  parameters = list(
    surv_par = 0,      # Logit persistence probability
    surv_re = rep(0,max(d$grp)),
    detect_par = 0,
    detect_re = rep(0,max(d$grp)),
    taggingRate_par = rep(0,max(d$grp)),      # Logit tagging probability
    Nsuper = 6,   # Log expected carcass births
    Nsuper_re = rep(0,max(d$grp)+1),   # Log expected carcass births
    Nsuper_re_sig = 0,
    pent_sig = rep(0,1),
    # pent_sig = rep(0,max(d$grp)+1),
    pent = array(0, c(max(d$grp)+1,max(d$t_k))),
    detect_re_sig = 0,
    surv_re_sig = 0
  )
  #
  # # Define RTMB Model
  rtmb_model <- function(parms){

    RTMB::getAll(data, parms)

    nll <- 0

    max_l <- max(na.omit(r_l)) + 1

    # Compute probability matrix
    surv <- array(0,c(max(grp),2,2), dimnames = list("grp"=NULL,"init_state"=NULL,"state"=NULL))
    f_surv <- array(0,c(max(grp),2,2), dimnames = list("grp"=NULL,"init_state"=NULL,"state"=NULL))
    for(i in 1:max(grp)){
      surv[i,1,2] <- exp(surv_par + surv_re[i])
      surv[i,1,1] <- -surv[i,1,2]
      f_surv[i,,] <- surv[i,,]
      surv[i,,] <- as.matrix(Matrix::expm(surv[i,,]))
    }
    nll_surv_re <- RTMB::dnorm(surv_re, 0 , exp(surv_re_sig), log = TRUE)
    REPORT(surv)


    # Compute probability matrix
    p <- array(0,c(max(grp),2,2), dimnames = list("grp"=NULL,"init_state"=NULL,"state"=NULL))
    for(i in 1:max(grp)){
      p[i,1,2] <- exp(detect_par + detect_re[i])
      p[i,1,1] <- -p[i,1,2]
      p[i,,] <- as.matrix(Matrix::expm(p[i,,]))
    }
    nll_detect_re <- RTMB::dnorm(detect_re, 0 , exp(detect_re_sig), log = TRUE)
    REPORT(p)

    B <- exp(pent + Nsuper_re + Nsuper )
    Btotal <- sum(B)
    nll_pent <- RTMB::dnorm(pent, 0, exp(pent_sig), log = TRUE)
    nll_Nsuper <- RTMB::dnorm(Nsuper_re, 0, exp(Nsuper_re_sig), log = TRUE)

    REPORT(B)
    REPORT(pent)
    REPORT(Btotal)
    ADREPORT(Btotal)

    # Derived variables
    N <- array(0, c(max(grp)+1,2, max(t_k)),
               dimnames = list("grp"=NULL,"state"=NULL,"week"=NULL))      # Total available carcasses over time
    for(g in 1:(max(grp)+1)){
      if(g==1){
        tmp_mat <- t(as.matrix(Matrix::expm((f_surv[1,,] * B[2,1]  + f_surv[3,,] * B[4,1])/(B[2,1] + B[4,1]) * 0.5 )))
        N[g,,1] <- tmp_mat  %*% c(B[g,1],0)
        # Loop through time steps to compute derived variables
        for(t in 2:max(t_k)){
          tmp_mat <- t(as.matrix(Matrix::expm((f_surv[1,,] * B[2,t]  + f_surv[3,,] * B[4,t])/(B[2,t] + B[4,t]) * 0.5 )))
          N[g,,t] = tmp_mat %*% c(B[g,t],0) +
            t(surv[1,,]) %*% N[g,,t-1]
        }
      }else{
        N[g,,1] <- t(as.matrix(Matrix::expm(f_surv[g-1,,] * 0.5))) %*% c(B[g,1],0)
        # Loop through time steps to compute derived variables
        for(t in 2:max(t_k)){
          N[g,,t] = t(as.matrix(Matrix::expm(f_surv[g-1,,] * 0.5))) %*% c(B[g,t],0) +
            t(surv[g-1,,]) %*% N[g,,t-1]
        }
      }
    }
    # REPORT(N)
    print("test3")
    #
    nll_detect <- matrix(1,max(grp)+1,max(t_k))
    nll_tag <- matrix(0,max(grp)+1,max(t_k))
    TotalCarcasses_t <- array(0, c(max(grp)+1, max(t_k)), dimnames = list("grp"=NULL,"week"=NULL)) #Total
    CarcassesW_tags_t <- array(0, c(max(grp)+1, max(t_k)), dimnames = list("grp"=NULL,"week"=NULL)) #tags
    E_TotalCarcasses_t <- array(0, c(max(grp)+1, max(t_k)), dimnames = list("grp"=NULL,"week"=NULL)) #predicted
    E_TaggedCarcasses_t <- array(0, c(max(grp)+1, max(t_k)), dimnames = list("grp"=NULL,"week"=NULL)) #predicted
    for(g in 1:(max(grp)+1)){
      for (t in 1:max(t_k)) {
        # Detection Process (Binomial likelihood)
        CarcassesW_tags_t[g, t] <- sum(na.omit(n[grp == (g-1) & t_k == t & t_l == 1 & tag == 1]))  # Tagged carcasses at time t
        TotalCarcasses_t[g, t] <- sum(na.omit(n[grp == (g-1) & t_k == t & t_l == 1]))  # Total detected carcasses at time t
        #Carcass detection probability likelihood
        if(g==1){
          E_TotalCarcasses_t[g,t] <- t(c(1,0)) %*% ((t(p[1,,]) + t(p[3,,]))/2) %*% as.vector(N[g,,t])
        }else{
          E_TotalCarcasses_t[g,t] <- t(c(1,0)) %*% t(p[g-1,,]) %*% as.vector(N[g,,t])
        }
        nll_detect[g,t] <- RTMB::dpois(TotalCarcasses_t[g,t],
                                       E_TotalCarcasses_t[g,t] + 1,
                                       log = TRUE)
        if(g>1){
          #Carcass tagging rate likelihood
          E_TaggedCarcasses_t[g,t] <- E_TotalCarcasses_t[g, t] * plogis(taggingRate_par[g-1])

          nll_tag[g,t] <- RTMB::dpois(CarcassesW_tags_t[g,t],
                                      E_TaggedCarcasses_t[g,t] + 1,
                                      log = TRUE)
        }
      }

    }
    REPORT(TotalCarcasses_t)
    REPORT(E_TotalCarcasses_t)
    REPORT(E_TaggedCarcasses_t)
    REPORT(CarcassesW_tags_t)
    REPORT(nll_tag)
    REPORT(nll_detect)

    nll_CJS <- rep(0,length(t_k))
    for(i in 1:length(t_l)){
      m <- matrix(0,max_l,max_l)
      diag(m) <- 1
      last_k <- ifelse(is.na(r_k[i]),max(na.omit(r_k)),r_k[i])

      if(grp[i]>0){
        if(is.na(r_k[i])){#never recapped
          if((t_k[i]+1)<last_k){
            for(k in (t_k[i]+1):last_k){
              m <- m %*% surv[grp[i],,] %*% diag(p[grp[i],,max_l])
            }
          }
        }


        if(!is.na(r_k[i])){#Recaptured
          if(r_k[i]>0){
            for(k in (t_k[i]+1):r_k[i]){
              if(k<last_k){
                m <- m %*% surv[grp[i],,] %*% diag(p[grp[i],,max_l]) #non-detection
              }else{
                m <- m %*% surv[grp[i],,] %*% diag(p[grp[i],,r_l[i]]) #recap
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
      }else{
        nll_CJS[i] <- 1
      }
    }

    # pred_n <- nll_CJS *
    nll_CJS_total <- sum(log(nll_CJS[tag==TRUE])*n[tag==TRUE])
    REPORT(nll_CJS_total)


    return(-(nll_CJS_total +
               sum(nll_detect_re) +
               sum(nll_surv_re) +
               sum(nll_Nsuper) +
               sum(nll_detect) +
               sum(nll_tag) +
               sum(nll_pent))
    )
    # return(0)
  }
  #
  obj <- RTMB::MakeADFun(rtmb_model,
                         parameters,
                         random = c("pent","surv_re","detect_re","Nsuper_re"),
                         map = list(
                           Nsuper_re = as.factor(rep(NA, length(parameters$Nsuper_re)))
                           ,Nsuper_re_sig = as.factor(NA)
                           # surv_re_sig = as.factor(NA)
                                    # ,detect_re = as.factor(rep(NA, length(parameters$detect_re)))
                                    # ,detect_re_sig = as.factor(NA)
                         ),
                         silent = FALSE)

  # Optimize Model
  opt <- TMBhelper::fit_tmb(obj)
  rep <- obj$report()
  sdr <- sdreport(obj)

  lk_up <- d %>%
    group_by(grp,grp_id) %>%
    mutate(grp = grp + 1) %>%
    summarise(n = sum(n))


  data.frame(reshape2::melt(rep$E_TotalCarcasses_t),
             obs = reshape2::melt(rep$TotalCarcasses_t)$value) %>%
    left_join(lk_up,"grp") %>%
    mutate(year = as.factor(substr(grp_id,1,4)),
           fish_type = substr(grp_id,6,length(grp_id))) %>%
    ggplot(aes(x = week, y = obs, color = year)) +
    geom_point() +
    geom_line(aes(x =  week, y = value, color = year)) +
    facet_wrap(~fish_type, ncol = 2, scales = "free") +
    labs(title = "Total carcasses observed")

  data.frame(reshape2::melt(rep$E_TaggedCarcasses_t),
             obs = reshape2::melt(rep$CarcassesW_tags_t)$value) %>%
    left_join(lk_up,"grp") %>%
    mutate(year = as.factor(substr(grp_id,1,4)),
           fish_type = substr(grp_id,6,length(grp_id))) %>%
    ggplot(aes(x = week, y = obs, color = year)) +
    geom_point() +
    geom_line(aes(x =  week, y = value, color = year)) +
    facet_wrap(~fish_type, ncol = 2, scales = "free") +
    labs( title = "Total number of tags")
  #
  data.frame(reshape2::melt(rep$p)) %>%
    filter(init_state <= state) %>%
    mutate(init_state = as.factor(init_state),
           state = factor(state)) %>%
    mutate(grp = grp + 1) %>%
    left_join(lk_up,"grp") %>%
    ggplot(aes(x = state, y = init_state, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    facet_wrap(~grp_id, ncol = 2, scales = "free") +
    geom_text(aes(label = round(value, 2)),
              color = "black",
              size = 5) +
    theme(text = element_text(size = 14)) +
    labs(title = "Detection probability")


  data.frame(reshape2::melt(rep$surv)) %>%
    filter(init_state <= state) %>%
    mutate(init_state = as.factor(init_state),
           state = factor(state)) %>%
    mutate(grp = grp + 1) %>%
    left_join(lk_up,"grp") %>%
    ggplot(aes(x = state, y = init_state, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    facet_wrap(~grp_id, ncol = 2, scales = "free") +
    geom_text(aes(label = round(value, 2)),
              color = "black",
              size = 5) +
    theme(text = element_text(size = 14)) +
    labs(title = "Survival")


  print(rep$B)
  print(year)
  print(rep$Btotal)
